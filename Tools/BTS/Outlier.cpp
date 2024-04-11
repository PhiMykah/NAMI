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
        vector csim = vector(CSimMSD(matrix, n_atoms));
        vector::iterator it = min_element(csim.begin(), csim.end());
        int index = std::distance(csim.begin(), it);
        return index;
    }
    int N = matrix.N;
    Matrix sq_data_total = matrix.pow(2);
    vector c_sum_total = matrix.Sum(COL);
    vector sq_sum_total = sq_data_total.Sum(COL);
    int index = N + 1;
    int max_dissim = INT_MAX;

    for (int row = 0; row < N; row++){
        float value = ExtendedComparison(c_sum_total - matrix[row], metric, N-1, n_atoms);
        if (value < max_dissim) {
            max_dissim = value;
            index = row;   
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
    
    int N = matrix.N;
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
    
    int N = matrix.N;
    int cutoff = n_trimmed;

    if (criterion == Criterion::SIM_TO_MEDOID) {
        int medoid_index = CalculateMedoid(matrix, metric, n_atoms);
        vector medoid = matrix[medoid_index];
        // Remove the values from the medoid index of matrix
        matrix.erase(medoid_index);
    
        vector values;
        for (int row = 0; row < N; row++){
            values.push_back(
                ExtendedComparison(matrix[row], medoid, metric, n_atoms) // data_type = full?
            );
        }

        // Sort the list
        std::vector<int> highest_indices;
        for (long unsigned int i = 0; i < values.size(); i++)
        {
            highest_indices.push_back(i);
        }

        quicksort(values, highest_indices, 0, values.size()-1);

        // Collect the indices of the last cutoff elements of the sorted list
        highest_indices.erase(highest_indices.end()-cutoff, highest_indices.end()); 

        // Remove the values at those indices of the original matrix
        Matrix newMatrix(matrix.GetArray());

        newMatrix.erase(highest_indices);

        return newMatrix;

    } else {
        vector c_sum = matrix.Sum(COL);
        vector sq_sum_total = matrix.pow(2).Sum(COL);
        vector comp_sims;
        vector values;
        for (int row = 0; row < N; row++){
            vector c = c_sum - matrix[row];
            vector sq = sq_sum_total - (pow(matrix[row],2));
            values.push_back(
                ExtendedComparison(c, sq, metric, N-1, n_atoms)
            );
        }
        // Sort the list
        std::vector<int> lowest_indices;
        for (long unsigned int i = 0; i < values.size(); i++)
        {
            lowest_indices.push_back(i);
        }

        quicksort(values, lowest_indices, 0, values.size()-1);

        // Collect the indices of the first cutoff elements of the sorted list
        lowest_indices.erase(lowest_indices.begin()+cutoff, lowest_indices.end()); 
        
        // Remove the values at those indices of the original matrix
        Matrix newMatrix(matrix.GetArray());

        newMatrix.erase(lowest_indices);

        return newMatrix;
    }

}
