import numpy as np
from mdance import extended_comparison, calculate_comp_sim
from mdance.tools.bts import calculate_medoid, calculate_outlier
from mdance.tools.bts import trim_outliers, diversity_selection
from mdance.tools.bts import align_traj, equil_align

FILE = "examples/backbone.npy"
metrics = ['MSD', 'BUB', 'Fai', 'Gle', 'Ja', 'JT', 'RT', 'RR', 'SM', 'SS1', 'SS2']

def output_results(title : str, results : list):
    print(f"{title}:")
    for i in range(len(results)):
        pad = f" " if len(metrics[i]) == 2 else ""
        print(f"Metric {pad}\'{metrics[i]}\': {results[i]}")
    print("\n")

def main():
    matrix = np.load(FILE)
    int_trimmed = 3
    float_trimmed = .5
    N_atoms = 10
    percentage = 75

    ec_vals = [extended_comparison(matrix, metric=metric, N_atoms=N_atoms) for metric in metrics]

    ccs_vals = [calculate_comp_sim(matrix, metric=metric, N_atoms=N_atoms) for metric in metrics]

    cm_vals = [calculate_medoid(matrix, metric=metric, N_atoms=N_atoms) for metric in metrics]

    co_vals = [calculate_outlier(matrix, metric=metric, N_atoms=N_atoms) for metric in metrics]

    toi_vals = [trim_outliers(matrix, int_trimmed, metric, N_atoms=N_atoms) for metric in metrics]

    tof_vals = [trim_outliers(matrix, float_trimmed, metric, N_atoms=N_atoms) for metric in metrics]

    tofm_vals = [trim_outliers(matrix, float_trimmed, metric, N_atoms=N_atoms, criterion='sim_to_medoid') for metric in metrics]

    # dsm_vals = [diversity_selection(matrix, percentage, metric, start='medoid', N_atoms=N_atoms) for metric in metrics]

    # dso_vals = [diversity_selection(matrix, percentage, metric, start='outlier', N_atoms=N_atoms) for metric in metrics]

    # dsr_vals = [diversity_selection(matrix, percentage, metric, start='random', N_atoms=N_atoms) for metric in metrics]
    
    # dsl_vals = [diversity_selection(matrix, percentage, metric, start=[0, 2, 4], N_atoms=N_atoms) for metric in metrics]  # type: ignore

    # atu_vals = [align_traj(matrix, N_atoms, align_method='uni')]

    # atk_vals = [align_traj(matrix, N_atoms, align_method='kron')]

    tests : list[tuple] = [
        ("Extended Comparison", ec_vals),
        ("Calculate Complementary Similarity", ccs_vals),
        ("Calculate Medoid", cm_vals),
        ("Calculate Outlier", co_vals),
        ("Trim Outliers Integer", toi_vals),
        ("Trim Outliers Float", tof_vals),
        ("Trim Outliers Float (MEDOID)", tof_vals),
        # ("Diversity Selection + NewIndex (Medoid)", dsm_vals),
        # ("Diversity Selection + NewIndex (Outlier)", dso_vals),
        # ("Diversity Selection + NewIndex (Random)", dsr_vals),
        # ("Diversity Selection + NewIndex (0, 2 , 4)", dsl_vals),
        # ("Align Trajectory Uniform", atu_vals),
        # ("Align Trajectory Kronecker", atk_vals),
    ]

    for test in tests:
        output_results(test[0], test[1])


if __name__ == "__main__":
    main()