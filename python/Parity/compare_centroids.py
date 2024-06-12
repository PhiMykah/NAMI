import numpy as np
import os

def get_dif(file1 : str, file2 : str, skip1 = 0, skip2 = 0):
    a = np.loadtxt(file1, delimiter=',', skiprows=skip1)

    b = np.loadtxt(file2, delimiter=',', skiprows=skip2)

    dif = a - b

    return (a.min(), a.max(), b.min(), b.max(), dif.min(), dif.max())

def main():
    PATH1 = f"{os.environ['NAMI']}/output/ala10/NAMI"
    PATH2 = f"{os.environ['NAMI']}/output/ala10/MDANCE"

    initiators = ['random_centroids', 'comp_sim_centroids', 'div_sel_centroids']
    confirmations = ['alpha', 'hairpin', 'left', 'pp2']

    tests = []
    for conf in confirmations:
        for init in initiators:
            tests.append(f"{init}_{conf}")

    nami_files = [f"{PATH1}/{test}.csv" for test in tests]
    mdance_files = [f"{PATH2}/{test}.csv" for test in tests]

    for n,m,t in zip(nami_files, mdance_files, tests):
        n_min, n_max, m_min, m_max, d_min, d_max = get_dif(n,m, skip1=1, skip2=4)
        print(f"{t}\n" + "-" * len(t))
        print(f"N - Min:\t{n_min: 8.4f}\tMax:\t{n_max: 8.4f}")
        print(f"M - Min:\t{m_min: 8.4f}\tMax:\t{m_max: 8.4f}")
        print(f"D - Min:\t{d_min: 8.4f}\tMax:\t{d_max: 8.4f}")
        print("\n\n")

if __name__ == "__main__":
    main()
