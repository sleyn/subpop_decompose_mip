import os
import numpy as np
import pandas as pd
from collections import defaultdict
import time

# set tests
test_dir = os.path.join('Test', 'Test_huber_high_lambda')
os.mkdir(test_dir)


def p_clon_init(prob, n_variants, n_subpop):
    print('Initialize Presence of variants matrix')
    test_prob = np.random.choice([0, 1], (n_variants, n_subpop), p=[1 - prob, prob])
    # If there is WT - add variant to it
    test_prob[np.random.randint(n_variants), np.sum(test_prob, axis=0) == 0] = 1
    return test_prob


def f_clon_init(dirichlet_papam, n_subpop, n_samples):
    print('Initialize Sub-population frequency matrix matrix')
    f_clon = np.zeros((n_subpop, n_samples))
    for i in range(n_samples):
        f_clon[:, i] = np.random.dirichlet(np.ones(n_subpop) * dirichlet_papam, size=1)

    return f_clon


def collapse_identical(p_clon):
    # Sum identical subpopulations
    drop_subpop = set()

    for subpop in range(p_clon.shape[1]):
        for subpop1 in range(subpop + 1, p_clon.shape[1]):
            if np.min(p_clon[:, subpop] == p_clon[:, subpop1]):
                drop_subpop.add(subpop1)

    p_clon = np.delete(p_clon, list(drop_subpop), axis=1)
    return p_clon

number_of_tests = 30
results_df = pd.DataFrame({
    'N': np.arange(0, number_of_tests, 1),     # Number of the test
    'Variants': np.zeros(number_of_tests),
    'Samples': np.zeros(number_of_tests),
    'Subpopulations': np.zeros(number_of_tests),
    'TP': np.zeros(number_of_tests),
    'FP': np.zeros(number_of_tests),
    'FN': np.zeros(number_of_tests),
    'Precision': np.zeros(number_of_tests),
    'Recall': np.zeros(number_of_tests),
    'F1_score': np.zeros(number_of_tests),
    'MSE_SQRT': np.zeros(number_of_tests),            # Mean Square Error for Subpopulation matrix
    'Exec_Time': np.zeros(number_of_tests)
})

for i in range(number_of_tests):
    os.mkdir(os.path.join(test_dir, 'test' + str(i)))
    n_variants = np.random.randint(2, 100, 1)[0]
    n_samples = np.random.randint(1, 10, 1)[0]
    n_subpop = n_samples * 2
    test_p_clon = p_clon_init(0.1, n_variants, n_subpop)
    test_p_clon = collapse_identical(test_p_clon)
    n_subpop = test_p_clon.shape[1]
    test_f_clon = f_clon_init(0.1, n_subpop, n_samples)
    test_f_var = np.dot(test_p_clon, test_f_clon)
    test_f_var_df = pd.DataFrame(test_f_var)
    test_f_var_df['Chrom'] = 'X'
    test_f_var_df['Pos'] = 1
    test_f_var_df['Ref'] = 'A'
    test_f_var_df['Alt'] = 'T'

    print(f'Test run:\n\tNumber of varinats: {n_variants}\n\tNumber of subpopulations: {n_subpop}')
    test_f_var_df.to_csv(os.path.join(test_dir, 'test' + str(i), 'test_variant_table.tsv'), header=True,
                                    index=False, sep='\t')
    pd.DataFrame(test_p_clon).to_csv(os.path.join(test_dir, 'test' + str(i), 'test_presence_table.tsv'), header=False,
                                    index=False, sep='\t')
    pd.DataFrame(test_f_clon).to_csv(os.path.join(test_dir, 'test' + str(i), 'test_subpop_table.tsv'), header=False,
                                    index=False, sep='\t')

    print('Run ClonalDecomposition.')
    start_time = time.time()
    os.system(f'python decompose.py -t {os.path.join(test_dir, "test" + str(i), "test_variant_table.tsv")} '
              f'-o {os.path.join(test_dir, "test" + str(i), "out")} '
              f'&> {os.path.join(test_dir, "test" + str(i), "log.txt")}')
    exec_time = time.time() - start_time
    print(f'Runtime: {exec_time}')

    # To calculate:
    # 1. how many subpopulation profiles are correct (exclude WT)
    # 2. how many profiles are missed
    # 3. how many profiles are wrong
    # 4. distance between initial subpop matrix and reconstructed
    # 5. distance between initial variant matrix and reconstructed

    # Read reconstructed profile
    print('Calculate test statistics.')
    reconstructed_p = np.loadtxt(os.path.join(test_dir, "test" + str(i), "out", "presence_sub_collapsed.tsv"),
                                 dtype='i',
                                 delimiter='\t')

    # Remove WT from reconstructed profile
    wt_index = np.where(np.sum(reconstructed_p, axis=0) == 0)
    reconstructed_p = reconstructed_p[:, np.sum(reconstructed_p, axis=0) > 0]

    reconstructed_p_df = pd.DataFrame(reconstructed_p.T)
    reconstructed_p_df['profile'] = reconstructed_p_df.apply(
        lambda x: ''.join([str(presence) for presence in x.tolist()]),
        axis=1)

    # Read reconstructed subpopulation matrix
    reconstructed_s = np.loadtxt(os.path.join(test_dir, "test" + str(i), "out", "freq_sub_collapsed.tsv"),
                                 delimiter='\t')

    # Remove WT
    reconstructed_s = np.delete(reconstructed_s, wt_index, 0)

    reconstructed_df = pd.concat([reconstructed_p_df['profile'], pd.DataFrame(reconstructed_s)], axis=1)

    # Create initial data frame
    initial_p_df = pd.DataFrame(test_p_clon.T)
    initial_p_df['profile'] = initial_p_df.apply(
        lambda x: ''.join([str(presence) for presence in x.tolist()]),
        axis=1)

    initial_df = pd.concat([initial_p_df['profile'], pd.DataFrame(test_f_clon)], axis=1)

    comparison_df = initial_df.merge(reconstructed_df, how='outer', on='profile')
    initial_col_index = 1
    reconstructed_col_index = initial_df.shape[1]

    def checkConsistency(x, init_ind, recon_ind):
        if x.isnull()[init_ind]:
            return 'FP'
        if x.isnull()[recon_ind]:
            return 'FN'
        return 'TP'


    comparison_df['Rec_type'] = comparison_df.apply(
        lambda x: checkConsistency(x, initial_col_index, reconstructed_col_index),
        axis=1
    )

    counts = comparison_df['Rec_type'].value_counts()
    counts_dict = defaultdict(int, list(zip(counts.index.to_list(), counts.to_list())))
    tp = counts_dict['TP']
    fp = counts_dict['FP']
    fn = counts_dict['FN']

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    if precision > 0 or recall > 0:
        f1 = 2 * (precision * recall) / (precision + recall)
    else:
        f1 = 0

    comparison_df = comparison_df.fillna(0)
    init_s_combined = comparison_df.iloc[:, 1:(n_samples+1)].to_numpy()
    reconstructed_s_combined = comparison_df.iloc[:, (n_samples+1):comparison_df.shape[1]-1].to_numpy()
    mse = np.sqrt(1 / init_s_combined.size * np.sum(np.power(init_s_combined - reconstructed_s_combined, 2)))

    results_df.loc[i, 'N'] = i
    results_df.loc[i, 'Variants'] = n_variants
    results_df.loc[i, 'Samples'] = n_samples
    results_df.loc[i, 'Subpopulations'] = n_subpop
    results_df.loc[i, 'TP'] = tp
    results_df.loc[i, 'FP'] = fp
    results_df.loc[i, 'FN'] = fn
    results_df.loc[i, 'Precision'] = precision
    results_df.loc[i, 'Recall'] = recall
    results_df.loc[i, 'F1_score'] = f1
    results_df.loc[i, 'MSE_SQRT'] = mse
    results_df.loc[i, 'Exec_Time'] = exec_time
    print(results_df.loc[i, :])

    comparison_df.to_csv(os.path.join(test_dir,  "test" + str(i), 'Comparison_table.tsv'), index=False, sep='\t')


for col in ['N', 'Variants', 'Samples', 'Subpopulations', 'TP', 'FP', 'FN']:
    results_df[col] = round(results_df[col], 0)

for col in ['Precision', 'Recall', 'F1_score', 'MSE_SQRT', 'Exec_Time']:
    results_df[col] = round(results_df[col], 3)

results_df.to_csv(os.path.join(test_dir, 'Test_results.tsv'), index=False, sep='\t')
