#include "Outlier.h"

/*
Calculates the outliers of a dataset using the metrics in extended comparison.
Outliers are the least representative objects of a set.

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
    The index of the outlier in the dataset.
*/
int CalculateOutlier(Matrix matrix, Metric metric, int n_atoms){
    if (metric == Metric::MSD) {
        // Returns the indices where the minimum value occurs
        rvector csim = CSimMSD(matrix, n_atoms).as_row();
        uword min = csim.index_min();
        return min;
    }
    uword N = matrix.n_rows;
    Matrix sq_data_total = arma::pow(matrix,2);
    rvector c_sum_total = arma::sum(matrix,COL);
    rvector sq_sum_total = arma::sum(sq_data_total,COL);
    rvector c_diff;
    rvector sq_diff;
    uword index = N + 1;
    float min_dissim = INFINITY;
    float value;
    for (uword i = 0; i < N; i++){
        c_diff = c_sum_total - matrix.row(i);
        sq_diff = sq_sum_total - arma::pow(matrix.row(i), 2);
        value = ExtendedComparison(c_diff, sq_diff, metric, N-1, n_atoms);
        if (value < min_dissim) { // Might need to specify a significant difference between these two for assurance - φ
            min_dissim = value;
            index = i;   
        }
    }

    return index;

}


/*
Trims a desired percentage of outliers (most dissimilar) from the dataset 
by calculating largest complement similarity.

Parameters
----------
matrix : Matrix
    Input data matrix.
percent_trimmed : float
    The desired fraction of outliers to be removed.
metric : {'MSD', 'RR', 'JT', 'SM', etc}
    Metric used for extended comparisons. See `extended_comparison` for details.
N_atoms : int
    Number of atoms in the system.
criterion : {'comp_sim', 'sim_to_medoid'}, optional
    Criterion to use for data trimming. Defaults to 'comp_sim'.
    'comp_sim' removes the most dissimilar objects based on the complement similarity.
    'sim_to_medoid' removes the most dissimilar objects based on the similarity to the medoid.
    
Returns
-------
Matrix
    A matrix with desired fraction of outliers removed.

Notes
-----
If the criterion is 'comp_sim', the lowest indices are removed because they are the most outlier.
However, if the criterion is 'sim_to_medoid', the highest indices are removed because they are farthest from the medoid.
*/
Matrix TrimOutliers(
    Matrix matrix, float percent_trimmed, Metric metric, 
    int n_atoms, Criterion criterion){
    
    int N = matrix.n_rows;
    int cutoff = int(floor(N * percent_trimmed));
    
    return TrimOutliers(matrix, cutoff, metric, n_atoms, criterion);
}

/*
Trims a certain amount (most dissimilar) from the dataset 
by calculating largest complement similarity.

Parameters
----------
matrix : Matrix
    Input data matrix.
n_trimmed : int
    The desired number of outliers to be removed.
metric : {'MSD', 'RR', 'JT', 'SM', etc}x
    Metric used for extended comparisons. See `extended_comparison` for details.
N_atoms : int
    Number of atoms in the system.
criterion : {'comp_sim', 'sim_to_medoid'}, optional
    Criterion to use for data trimming. Defaults to 'comp_sim'.
    'comp_sim' removes the most dissimilar objects based on the complement similarity.
    'sim_to_medoid' removes the most dissimilar objects based on the similarity to the medoid.
    
Returns
-------
Matrix
    A matrix with desired fraction of outliers removed.

Notes
-----
If the criterion is 'comp_sim', the lowest indices are removed because they are the most outlier.
However, if the criterion is 'sim_to_medoid', the highest indices are removed because they are farthest from the medoid.
*/
Matrix TrimOutliers(
    Matrix matrix, int n_trimmed, Metric metric,
    int n_atoms, Criterion criterion){
    
    uword N = matrix.n_rows;
    int cutoff = n_trimmed;

    if (criterion == Criterion::SIM_TO_MEDOID) {
        int medoid_index = CalculateMedoid(matrix, metric, n_atoms);
        rvector medoid = matrix.row(medoid_index);
        // Remove the values from the medoid index of matrix
        matrix.shed_row(medoid_index);
    
        rvector values(N, arma::fill::zeros);
        for (uword i = 0; i < matrix.n_rows; i++){
            values(i) = ExtendedComparison(matrix.row(i), medoid, metric, n_atoms); // data_type = full?
        }

        // Sort the values
        index_vec sorted_indices = arma::sort_index(values);

        // Collect the indices of the last cutoff elements of the sorted list
        index_vec highest_indices = sorted_indices.subvec(sorted_indices.size()-cutoff, sorted_indices.size()-1);

        // Remove the values at those indices of the original matrix
        Matrix newMatrix(matrix);
    
        newMatrix.shed_rows(highest_indices);

        return newMatrix;

    } else {
        rvector c_sum = arma::sum(matrix,COL);
        rvector sq_sum_total = arma::sum(arma::pow(matrix,2),COL);
        rvector comp_sims;
        float result;
        rvector values (N, arma::fill::zeros);
        for (uword i = 0; i < N; i++){
            rvector c = c_sum - matrix.row(i);
            rvector sq = sq_sum_total - (arma::pow(matrix.row(i),2));
            result = ExtendedComparison(c, sq, metric, N-1, n_atoms);
            values(i) = result;
        }

        // Sort the indices of the values
        index_vec sorted_indices = arma::sort_index(values);

        // Collect the indices of the first cutoff elements of the sorted list
        index_vec lowest_indices = sorted_indices.subvec(0, cutoff);

        // Remove the values at those indices of the original matrix
        Matrix newMatrix(matrix);

        newMatrix.shed_rows(lowest_indices);

        return newMatrix;
    }

}
