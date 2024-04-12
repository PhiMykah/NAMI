#include "NewIndex.h"

/*
Function to get the new index to add to the selected indices

Parameters
----------
matrix : Matrix
    Input data matrix.
metric : {'MSD', 'RR', 'JT', 'SM', etc}
    Metric used for extended comparisons. See `extended_comparison` for details.
selected_condensed : vector
    Condensed sum of the selected fingerprints.
n : int
    Number of selected objects.
select_from_n : 1D array-like
    Indices of the objects to select from.
N_atoms : int, optional
    Number of atoms in the system. Defaults to 1.

Returns
-------
int
    index of the new fingerprint to add to the selected indices.
*/
int GetNewIndexN(Matrix matrix, Metric metric, vector select_condensed,
    int N, std::vector<int> select_from_n, int n_atoms)
{
    float sim_index;
    int n_total = N + 1;
    int min_value = INT_MIN;
    int index = matrix.N + 1;

    for (long unsigned int i = 0; i < select_from_n.size(); i++){
        // The extended comparison call may not be the right one here, check to see
        sim_index = ExtendedComparison(
            select_condensed + matrix[i], metric, n_total, n_atoms
            );
        if (sim_index > min_value) {
            min_value = sim_index;
            index = i;
        }
    }

    return index;
}

/*
Function to get the new index to add to the selected indices
with condensed sum sincluded

Parameters
----------
matrix : Matrix
    Input data matrix.
metric : {'MSD', 'RR', 'JT', 'SM', etc}
    Metric used for extended comparisons. See `extended_comparison` for details.
selected_condensed : vector
    Condensed sum of the selected fingerprints.
sq_selected_condensed : vector
    Condensed sum of the squared selected fingerprints. Defaults to None.
n : int
    Number of selected objects.
select_from_n : 1D array-like
    Indices of the objects to select from.

N_atoms : int, optional
    Number of atoms in the system. Defaults to 1.

Returns
-------
int
    index of the new fingerprint to add to the selected indices.
*/
int GetNewIndexN(Matrix matrix, Metric metric, vector select_condensed,
    vector sq_selected_condensed, int N, std::vector<int> select_from_n,
    int n_atoms)
{
    float sim_index;
    int n_total = N + 1;
    int min_value = INT_MIN;
    int index = matrix.N + 1;

    for (long unsigned int i = 0; i < select_from_n.size(); i++){
        // The extended comparison call may not be the right one here, check to see
        sim_index = ExtendedComparison(
            vector(select_condensed + matrix[i]), 
            vector(sq_selected_condensed + (pow(matrix[i],2))), 
            metric, n_total, n_atoms
            );

        if (sim_index > min_value) {
            min_value = sim_index;
            index = i;
        }
    }

    return index;
}