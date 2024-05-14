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
        vector csim = CSimMSD(matrix, n_atoms);
        uword max = csim.index_max();
        return max;
    }
    uword N = matrix.n_rows;
    Matrix sq_data_total = arma::pow(matrix,2);
    rvector c_sum_total = arma::sum(matrix,COL);
    rvector sq_sum_total = arma::sum(sq_data_total,COL);
    rvector diff;
    uword index = N + 1;
    int max_dissim = -1;

    for (uword i = 0; i < N; i++){
        diff = c_sum_total - matrix.row(i);
        float value = ExtendedComparison(diff, metric, N-1, n_atoms);
        if (value > max_dissim) { // Might need to specify a significant difference between these two for assurance - φ
            max_dissim = value;
            index = i;
        }
    }

    return (int)index;
}