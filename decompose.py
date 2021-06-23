from pulp import *
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='''
       Performs observed variant frequency matrix decomposition
       into variant presence binary matrix and subpopulation frequency
       matrces using mixed integer programming.
''')

parser.add_argument('-t', '--var_table', type=str,
                    help='Variant table. Columns: Chrom, Pos, Ref, Alt and sample columns.')
parser.add_argument('-p', '--perc2prop', action='store_true',
                        help='If frequencies are in percents bring them '
                             'to proportions. Default: False')
parser.add_argument('-o', '--out_dir', type=str, default='output',
                        help='Output directory. Default: output')
parser.add_argument('--timelimit', type=int, default=180,
                        help='Time limit for solving the problem. Default: 180')
args = parser.parse_args()

try:
    os.makedirs(args.out_dir)
except OSError:
    print(f'Can\'t create output directory {args.out_dir}')
    exit(1)

# Read observed variants file (o_matrix)
# read variant frequency table
freq_var = pd.read_csv(args.var_table, sep='\t')
freq_var_columns = freq_var.columns.to_list()
sample_columns = [column for column in freq_var_columns if column not in ['Chrom', 'Pos', 'Ref', 'Alt']]

# convert percents to proportions
if args.perc2prop:
    freq_var[sample_columns] = freq_var[sample_columns] / 100

freq_var_matrix = freq_var[sample_columns].to_numpy()

# Collect dimensions
n_variants = freq_var_matrix.shape[0]
n_subpop = freq_var_matrix.shape[0] * 2
n_samples = freq_var_matrix.shape[1]

o_matrix = freq_var_matrix.tolist()

# make solver for one sample
SUBPOP = range(n_subpop)
VARIANTS = range(n_variants)
SAMPLES = range(n_samples)

# Solver Gurobi with 20min time limit
solver = GUROBI_CMD(options=[('timeLimit', args.timelimit), ('msg', True)])

prob = LpProblem('Sample_clonal_decomposition', LpMinimize)
p_matrix = LpVariable.dicts('presence', (VARIANTS, SUBPOP), cat='Binary')
s_matrix = LpVariable.dicts('subpop', (SUBPOP, SAMPLES), lowBound=0, upBound=1)
z_matrix = LpVariable.dicts('multiplex', (VARIANTS, SUBPOP, SAMPLES), lowBound=0, upBound=1)

# Linearize multiplication of binary matrix P and subpopulation matrix P
# z1 + z2 + z3 = o11
# for i in 1:3:
# 	zi <= p1i
# 	zi >= si1 + p1i - 1
#   zi >= 0
#   zi <= si1
for m in SAMPLES:
    print(f'Sample {m}')
    for v in VARIANTS:
        prob += lpSum([z_matrix[v][s][m] for s in SUBPOP]) == o_matrix[v][m]
        for s in SUBPOP:
            prob += z_matrix[v][s][m] - p_matrix[v][s] <= 0
            prob += z_matrix[v][s][m] - s_matrix[s][m] - p_matrix[v][s] >= -1
            prob += z_matrix[v][s][m] >= 0
            prob += z_matrix[v][s][m] - s_matrix[s][m] <= 0

    prob += lpSum([s_matrix[s][m] for s in SUBPOP]) <= 1

# Solve the problem
print('Solve a problem')
prob.solve(solver)
print("Status:", LpStatus[prob.status])

# Output subpopulations matrix
with open(os.path.join(args.out_dir, 'freq_sub.tsv'), 'w') as out_s_file:
    for s in SUBPOP:
        # print(f'{s}: {value(s_matrix[s])}')
        out_s_file.write('\t'.join([str(value(s_matrix[s][m])) for m in SAMPLES]) + '\n')

# Output presence matrix
with open(os.path.join(args.out_dir, 'presence_sub.tsv'), 'w') as out_p_file:
    for v in VARIANTS:
        out_p_file.write('\t'.join([str(value(p_matrix[v][s])) for s in SUBPOP]) + '\n')

# Output reconstructd matrix
reconstructed = np.ones((n_variants, n_samples))
np_s = np.ones((n_subpop, n_samples))
np_p = np.ones((n_variants, n_subpop))

# convert p_matrix to numpy array
for v in VARIANTS:
    for s in SUBPOP:
        np_p[v, s] = value(p_matrix[v][s])

# convert s_matrix to numpy array
for s in SUBPOP:
    for m in SAMPLES:
        np_s[s, m] = value(s_matrix[s][m])
        

# Sum identical subpopulations
drop_subpop = set()

for subpop in range(np_p.shape[1]):
    # Additionally check if subpopulation has any frequency at least in one sample
    if np.sum(np_s[subpop, :]) == 0:
        drop_subpop.add(subpop)
    for subpop1 in range(subpop + 1, np_p.shape[1]):
        if np.min(np_p[:, subpop] == np_p[:, subpop1]):
            np_s[subpop, :] += np_s[subpop1, :]
            np_s[subpop1, :] = 0
            drop_subpop.add(subpop1)

np_p = np.delete(np_p, list(drop_subpop), axis=1)
np_s = np.delete(np_s, list(drop_subpop), axis=0)

# Add WT (all variants are 0) if not present and add up to 100%
wt_index = -1

for subpop in range(np_p.shape[1]):
    # Check if it is wild-type
    if np.min(np_p[:, subpop] == np.zeros((np_p.shape[0], 1))):
        wt_index = subpop

# Define RELU function
def relu(x):
    return np.maximum(0, x)

# sum up WT to 100%
if wt_index >= 0:
    np_s[wt_index, :] += np.array(list(map(relu, 1 - np.sum(np_s, axis=0))))
else:
    np_p = np.concatenate((np_p, np.zeros((np_p.shape[0], 1))), axis=1)
    np_s = np.concatenate((np_s, np.zeros((1, np_s.shape[1]))), axis=0)
    np_s[np_s.shape[0] - 1, :] += np.array(list(map(relu, 1 - np.sum(np_s, axis=0))))

# Print collapsed tables
pd.DataFrame(np_p).to_csv(os.path.join(args.out_dir, 'presence_sub_collapsed.tsv'), header=False, index=False, sep='\t')
pd.DataFrame(np_s).to_csv(os.path.join(args.out_dir, 'freq_sub_collapsed.tsv'), header=False, index=False, sep='\t')
pd.DataFrame(np_p.dot(np_s)).to_csv(os.path.join(args.out_dir, 'reconstructed_collapsed.tsv'),
                                                            header=False, index=False, sep='\t')

# Run parent-child determination
subpops_all = pd.DataFrame(np_p)
freqs_all = pd.DataFrame(np_s)

# Save number of subpopulations to the variable
n_subpop_final = len(subpops_all.columns)
sample_names = ['S' + str(i+1) for i in range(n_subpop_final)]
subpops_all = subpops_all.rename(columns=dict(zip(subpops_all.columns, sample_names)))
freqs_all['Sample'] = sample_names

def findParents(subpops, n_sp):
    # Technical table to write overall distance between two subpopulations and maximum observed distance
    # If vectorized Sub1 - Sub2 will result in maximum 0 for all elements, then Sub1 is parent of Sub2.
    distance_tbl = pd.DataFrame({
        'Parent': ['-'] * (n_sp ** 2 + 1),
        'Descendant': ['-'] * (n_sp ** 2 + 1),
        'distance': [0] * (n_sp ** 2 + 1),
        'maximum': [0] * (n_sp ** 2 + 1)
    })

    # Populate distance table
    for i in range(n_sp):
        for j in range(n_sp):
            if i == j:
                distance_tbl.iloc[n_sp * i + j, :] = [subpops.columns[i], subpops.columns[j], 0, 0]
                if np.sum(subpops.iloc[:, i].to_numpy()) == 0:
                    distance_tbl.iloc[-1, :] = [subpops.columns[i], subpops.columns[i], 1, 0]

            par = subpops.iloc[:, i].to_numpy()
            des = subpops.iloc[:, j].to_numpy()
            maximum = np.max(par - des)
            distance = np.sum(np.abs(par - des))

            distance_tbl.iloc[n_sp * i + j, :] = [subpops.columns[i], subpops.columns[j], distance, maximum]

    # Select parent-descendent pairs
    distance_tbl_zero = distance_tbl[(distance_tbl['maximum'] == 0) & (distance_tbl['distance'] > 0)]
    distance_tbl_zero = distance_tbl_zero.loc[distance_tbl_zero.groupby('Descendant')['distance'].idxmin()]

    return distance_tbl_zero


# Output table
distance_tbl_out = pd.DataFrame({
        'Parent': ['-'] * n_subpop_final,
        'Descendant': sample_names,
        'distance': [0] * n_subpop_final,
        'maximum': [0] * n_subpop_final
    })

# Samples for the cullrent round of parent evaluation
samples_to_analysis = set()

for i in range(n_samples):
    samples_to_analysis.update(
        freqs_all['Sample'][freqs_all.iloc[:, i] > 0].to_list()
    )
    subpops_to_analysis = subpops_all[samples_to_analysis]
    tbl_out_round = findParents(subpops_to_analysis, len(samples_to_analysis))

    for relation in tbl_out_round.iterrows():
        subpop_current = relation[1]['Descendant']
        # Write only if Parent not found yet
        if distance_tbl_out[distance_tbl_out['Descendant'] == subpop_current]['Parent'].item() == '-':
            distance_tbl_out.loc[distance_tbl_out['Descendant'] == subpop_current, 'Parent'] = relation[1]['Parent']
            distance_tbl_out.loc[distance_tbl_out['Descendant'] == subpop_current, 'distance'] = relation[1]['distance']
            distance_tbl_out.loc[distance_tbl_out['Descendant'] == subpop_current, 'maximum'] = relation[1]['maximum']


# Sort output
distance_tbl_out['number'] = [int(sample[1:]) for sample in distance_tbl_out['Descendant'].to_list()]
distance_tbl_out = distance_tbl_out.sort_values(by=['number']).drop(columns=['number'])

# Write output
distance_tbl_out.to_csv(os.path.join(args.out_dir, 'parent_collapsed.tsv'), sep='\t', index=False)
