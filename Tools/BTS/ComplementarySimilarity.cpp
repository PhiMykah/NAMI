#include "ComplementarySimilarity.h"

/*
Complementary similarity is calculating the similarity of a set
without one object or observation using metrics in the extended comparison.
The greater the complementary similarity, the more representative the object is.

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
Matrix
    Matrix of complementary similarities for each object.
*/
Matrix CalculateCompSim(Matrix matrix, Metric metric, int n_atoms){
    if (metric == Metric::MSD) {
        return CSimMSD(matrix, n_atoms);
    }

    int N = matrix.N;

    Matrix sq_data_total = matrix.pow(2);

    vector c_sum_total = matrix.Sum(COL);

    vector sq_sum_total = sq_data_total.Sum(COL); 

    vector values;
    for (int row = 0; row < N; row++){
        values.push_back(
            ExtendedComparison(c_sum_total - matrix[row], metric, N-1, n_atoms)
        );
    }

    return Matrix(vec2D {values});
}


// Simplified complementary similarity calculation if metric is MSD
Matrix CSimMSD(Matrix matrix, int n_atoms){
    int N = matrix.N;
    Matrix sq_data = matrix.pow(2);

    vector c_sum = matrix.Sum(COL);
    vector sq_sum = sq_data.Sum(COL);

    Matrix comp_csum = c_sum - matrix;
    Matrix comp_sqsum = sq_sum - sq_data;

    Matrix total = vec2D{(2 * ((N-1) * comp_sqsum - comp_csum.pow(2))).Sum(ROW)};

    Matrix comp_msd = total / std::pow(N-1, 2);

    Matrix norm_msd = comp_msd / n_atoms;

    return norm_msd;
}