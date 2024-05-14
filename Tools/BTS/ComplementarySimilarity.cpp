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
vector CalculateCompSim(Matrix matrix, Metric metric, int n_atoms){
    if ((metric == Metric::MSD) && (n_atoms == 1)) {
        fprintf(stderr, "n_atoms is being specified as 1. Please change if n_atoms is not 1.\n");
    } 
    if (metric == Metric::MSD) {
        return CSimMSD(matrix, n_atoms);
    }

    uword N = matrix.n_rows;

    Matrix sq_data_total = arma::pow(matrix,2);

    rvector c_sum_total = arma::sum(matrix,COL);

    rvector sq_sum_total = arma::sum(sq_data_total,COL); 
    rvector diff; 
    rvector values(N);
    for (uword i = 0; i < N; i++){
        diff = c_sum_total - matrix.row(i);
        values(i) =
            ExtendedComparison(diff, metric, N-1, n_atoms);
    }

    return values.t();
}


// Simplified complementary similarity calculation if metric is MSD
vector CSimMSD(Matrix matrix, int n_atoms){
    int N = matrix.n_rows;
    Matrix sq_data = arma::pow(matrix,2);

    rvector c_sum = arma::sum(matrix,COL);
    rvector sq_sum = arma::sum(sq_data,COL);

    Matrix comp_csum(matrix.n_rows, matrix.n_cols);
    
    for (uword i = 0; i < matrix.n_rows; i++) {
        comp_csum.row(i) = c_sum - matrix.row(i);
    }

    Matrix comp_sqsum(matrix.n_rows, matrix.n_cols);

    for (uword i = 0; i < matrix.n_rows; i++) {
        comp_sqsum.row(i) = sq_sum - sq_data.row(i);
    }

    vector total(arma::sum(2 * ((N-1) * comp_sqsum - arma::pow(comp_csum,2)),ROW));

    vector comp_msd = total / std::pow(N-1, 2);

    vector norm_msd = comp_msd / n_atoms;

    return norm_msd;
}