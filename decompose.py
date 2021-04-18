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

prob = LpProblem('Sample clonal decomposition', LpMinimize)
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
prob.solve()
print("Status:", LpStatus[prob.status])

# Output subpopulations matrix
with open(os.path.join(args.out_dir, 'freq_sub_collapsed.tsv'), 'w') as out_s_file:
    for s in SUBPOP:
        # print(f'{s}: {value(s_matrix[s])}')
        out_s_file.write('\t'.join([str(value(s_matrix[s][m])) for m in SAMPLES]) + '\n')

# Output presence matrix
with open(os.path.join(args.out_dir, 'presence_sub_collapsed.tsv'), 'w') as out_p_file:
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

pd.DataFrame(np_p.dot(np_s)).to_csv(os.path.join(args.out_dir, 'reconstructed_collapsed.tsv'),
                                                            header=False, index=False, sep='\t')