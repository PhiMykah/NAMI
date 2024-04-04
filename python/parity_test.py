import numpy as np
from mdance import extended_comparison, calculate_comp_sim

def csim(matrix, N_atoms = 1):
    N = len(matrix)
    sq_data = matrix ** 2
    c_sum = np.sum(matrix, axis=0)
    sq_sum = np.sum(sq_data, axis=0)
    comp_csum = c_sum - matrix
    comp_sqsum = sq_sum - sq_data
    comp_msd = np.sum(2 * ((N-1) * comp_sqsum - comp_csum ** 2), axis=1) / (N-1)**2
    norm_msd = comp_msd/N_atoms
    return norm_msd

if __name__ == "__main__":
    matrix = np.arange(1,26).reshape(5,5)
    
    N_atoms = 10
    metrics = ['MSD', 'BUB', 'Fai', 'Gle', 'Ja', 'JT', 'RT', 'RR', 'SM', 'SS1', 'SS2']

    ec_vals = [extended_comparison(matrix, metric=metric, N_atoms=N_atoms) for metric in metrics]

    ccs_vals = [calculate_comp_sim(matrix, metric=metric, N_atoms=N_atoms) for metric in metrics]

    norm_msd = csim(matrix, N_atoms)

    for i in range(len(ec_vals)):
        pad = f" " if len(metrics[i]) == 2 else ""
        print(f"Metric {pad}\'{metrics[i]}\': {ec_vals[i]}")

    for i in range(len(ccs_vals)):
        pad = f" " if len(metrics[i]) == 2 else ""
        print(f"Metric {pad}\'{metrics[i]}\':\n {ccs_vals[i]}")

