import numpy as np
import os
from pathlib import Path
from mdance import extended_comparison, calculate_comp_sim
from mdance.tools.bts import calculate_medoid, calculate_outlier
from mdance.tools.bts import trim_outliers, diversity_selection
from mdance.tools.bts import align_traj, equil_align
from mdance.modules.kmeansNANI.nani import KmeansNANI

FOLDER = f"{os.environ['NAMI']}/examples/ala10"
metric = 'MSD'
N_atoms = 109
N = 1000
n_clusters = 4
percent_trimmed = 0.8
percentage = 75
kmn_percentage = 10
n_iter = 10

kmeans_test_count = 3
initiator_titles = ["****************\nRandom Initiator\n****************\n",
             "******************\nComp Sim Initiator\n******************\n",
             "********************\nDiv Select Initiator\n********************\n"]
init_types = ['random', 'comp_sim', 'div_select']

out_files = [os.environ['NAMI'] + "/output/ala10/MDANCE/random_centroids_{}", os.environ['NAMI'] + "/output/ala10/MDANCE/comp_sim_centroids_{}", os.environ['NAMI'] + "/output/ala10/MDANCE/div_sel_centroids_{}"]

np.set_printoptions(threshold=20)

def main():
    files = sorted([f"{FOLDER}/{file}" for file in os.listdir(FOLDER)])
    fnames = sorted([f"{Path(file).stem}" for file in os.listdir(FOLDER)])
    
    for i in range(len(files)):
        matrix = np.load(files[i])
        print(f"\nConfirmation: {fnames[i]}\n---------------------\n")

        print("Extended Comparison\n{:.4f}".format(
            extended_comparison(matrix, metric=metric, N=N, N_atoms=N_atoms)
        ))

        print("Complementary Similarity\nResults:\n")
        print(calculate_comp_sim(matrix, metric=metric, N_atoms=N_atoms))

        print("Calculate Medoid\n{:d}".format(
            calculate_medoid(matrix, metric, N_atoms=N_atoms)
        ))

        print("Calculate Outlier\n{:d}".format(
            calculate_outlier(matrix, metric, N_atoms=N_atoms)
        ))

        print("Trim Outliers\nResults:\n")
        print(trim_outliers(matrix, percent_trimmed, metric, N_atoms, 'comp_sim'))

        print("Trim Outliers\nResults:\n")
        print(trim_outliers(matrix, percent_trimmed, metric, N_atoms, 'sim_to_medoid'))

        print("Diversity Selection + NewIndex (MEDOID)\n")
        ds = np.array(diversity_selection(matrix, percentage, metric, 'medoid', N_atoms), dtype=int)
        print(ds)
        
        print("------\nKMEANS\n------\n")
        for j in range(kmeans_test_count):
            print(initiator_titles[j])
            kmn = KmeansNANI(matrix, n_clusters, metric, N_atoms, init_types[j])
            if init_types[j] != 'random':
                initiators = kmn.initiate_kmeans()
            else:
                initiators = init_types[j]
                
            labels, centers, n_iter = kmn.kmeans_clustering(initiators)

            print(f"Number of iterations: {n_iter}")
            print("Labels:\n")
            print(labels)

            print("Centroids:\n")
            print(centers)
            
            output_file = out_files[j].format(fnames[i])

            Path(output_file).parent.mkdir(parents=True, exist_ok=True)

            header = f'Number of clusters: {n_clusters}, Number of iterations: {n_iter}\n\nCentroids\n'
            np.savetxt(f'{output_file}.csv', centers.T, delimiter=',', header=header)

if __name__ == "__main__":
    main()
