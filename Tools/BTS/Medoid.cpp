#include "Medoid.h"

/*
Calculates the medoid of a dataset using the metrics in extended comparison.
Medoid is the most representative object of a set.

Parameters
----------
matrix : Matrix
    Input data matrix.
metric : {'MSD', 'RR', 'JT', 'SM', etc}
    Metric used for extended comparisons. See `extended_comparison` for details.
N_atoms : int, optional
    Number of atoms in the system. Defaults to 1.

Returns
-------
int
    The index of the medoid in the dataset.
*/
int CalculateMedoid(Matrix matrix, Metric metric, int n_atoms){
    if (metric == Metric::MSD) {
        // Returns the indices where the maximum value occurs
        vector csim = vector(CSimMSD(matrix, n_atoms));
        vector::iterator it = max_element(csim.begin(), csim.end());
        int index = std::distance(csim.begin(), it);
        return index;
    }
    int N = matrix.N;
    Matrix sq_data_total = matrix.pow(2);
    vector c_sum_total = matrix.Sum(COL);
    vector sq_sum_total = sq_data_total.Sum(COL);
    int index = N + 1;
    int max_dissim = -1;

    for (int row = 0; row < N; row++){
        float value = ExtendedComparison(c_sum_total - matrix[row], metric, N-1, n_atoms);
        if (value > max_dissim) {
            max_dissim = value;
            index = row;
        }
    }

    return index;
}