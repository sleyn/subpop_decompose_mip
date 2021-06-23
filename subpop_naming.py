import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Creates names of subpopulations based on the variant table descriptions')
parser.add_argument('-p', '--p_table', help='Presence table. Output of the decompose.py script.')
parser.add_argument('-a', '--ann_table', help='Table with variant annotations.')
parser.add_argument('-o', '--out_table', help='Output table.')
parser.add_argument('-v', '--sep_var', default=' + ', help='Separation for enumeration of variants '
                                                           'in the subpopulation. Default: " + ".')
parser.add_argument('-s', '--sep_ann', default=':', help='Separation for enumeration of annotations of'
                                                         'a variant. Default: ":".')
parser.add_argument('-c', '--ann_cols', default='Gene Mutation', help='Columns from the variant annotation table '
                                                                      'that should be used for subpopulation '
                                                                      'annotaion. Variants should be in the same order'
                                                                      'and number as in the table that was decomposed '
                                                                      'by decompose.py script.'
                                                                      'Columns should be separated by space.')
args = parser.parse_args()

# Read presence table and name columns and variants
prs_table = pd.read_csv(args.p_table, sep='\t', header=None)
prs_table_colnames = ['S' + str(i+1) for i in range(prs_table.shape[1])]
prs_table = prs_table.rename(columns=dict(zip(prs_table.columns, prs_table_colnames)))
prs_n_v = prs_table.shape[0]

# Find WT to save it
sum_of_vars = prs_table.agg('sum')
WT_index = sum_of_vars.to_list().index(0) if 0 in sum_of_vars.to_list() else -1
if WT_index != -1:
    WT = sum_of_vars.index[sum_of_vars.to_list().index(0)]
else:
    WT = '-'

# Add names
prs_table['V_Name'] = ['V' + str(i+1) for i in range(prs_n_v)]



prs_table = pd.melt(prs_table, value_vars=[sp for sp in prs_table.columns if sp[0] == 'S'], id_vars=['V_Name']).\
    query('value > 0').\
    rename({'variabe': 'subpop'}).\
    drop(columns=['value']).\
    set_index('V_Name')

# Read annotation table and keep only selected columns
var_table = pd.read_csv(args.ann_table, sep='\t')
var_table = var_table[args.ann_cols.split(' ')]
var_n_v = var_table.shape[0]

if var_n_v == prs_n_v:
    var_table['V_Name'] = ['V' + str(i+1) for i in range(var_n_v)]
    var_table = var_table.set_index('V_Name')
else:
    print('Different number of variants in the annotation and decomposed presence tables. '
          'Variants should be in the same order and number as in the table that was decomposed '
          'by decompose.py script.')
    exit(1)

var_table['annotation'] = var_table.apply(lambda row: args.sep_ann.join(row.values.astype(str)), axis=1)
var_table = var_table.drop(columns=args.ann_cols.split(' '))

prs_table = prs_table.merge(var_table, how='left', on='V_Name').\
    groupby(['variable'])['annotation'].\
    apply(lambda x: args.sep_var.join(x)).\
    reset_index()

# Add WT back
if WT != '-':
    prs_table_WT = pd.DataFrame({'variable': [WT], 'annotation': ['WT']})
    prs_table = pd.concat([prs_table, prs_table_WT])

# Sort table
prs_table['number'] = [int(sample[1:]) for sample in prs_table['variable'].to_list()]
prs_table = prs_table.sort_values(by=['number']).drop(columns=['number'])

prs_table.\
    rename(columns={'variable': 'subpopulation'}).\
    to_csv(args.out_table, sep='\t', index=False)


