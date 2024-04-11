#include "DiversitySelection.h"

/* Selects a diverse subset of the data using the complementary similarity.

Parameters
----------
matrix : Matrix
    Input data matrix.
percentage : int
    Percentage of the data to select.
metric : {'MSD', 'RR', 'JT', 'SM', etc}
    Metric used for extended comparisons. See `extended_comparison` for details.
start : {'medoid', 'outlier', 'random'}, optional
    Seed of diversity selection. Defaults to 'medoid'.
N_atoms : int, optional
    Number of atoms in the system. Defaults to 1.

Returns
-------
list
    List of indices of the selected data.
*/
std::vector<int> DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    DiversitySeed start, int n_atoms)
{
    //int selected_n;
    int n_total = matrix.N;
    std::vector<int> total_indices;
    for (int i = 0; i < n_total; i++)
    {
        total_indices.push_back(i);
    }

    // if (start == DiversitySeed::OUTLIER){
    //     int seed = CalculateOutlier(matrix, metric, n_atoms);
    //     selected_n = seed;
    // } else if (start == DiversitySeed::RANDOM){
    //     int seed = rand() % n_total;
    //     selected_n = seed;
    // } else {
    //     int seed = CalculateMedoid(matrix, metric, n_atoms);
    //     selected_n = seed;
    // }

    //int N = 1;
     
    return total_indices; // Placeholder return
}

/* Selects a diverse subset of the data using the complementary similarity.

Parameters
----------
matrix : Matrix
    Input data matrix.
percentage : int
    Percentage of the data to select.
metric : {'MSD', 'RR', 'JT', 'SM', etc}
    Metric used for extended comparisons. See `extended_comparison` for details.
start : list
    Seed of diversity selection.
N_atoms : int, optional
    Number of atoms in the system. Defaults to 1.

Returns
-------
list
    List of indices of the selected data.
*/
std::vector<int> DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    std::vector<int> start, int n_atoms)
{
    std::vector<int> selected_n;
    int n_total = matrix.N;
    std::vector<int> total_indices;
    for (int i = 0; i < n_total; i++)
    {
        total_indices.push_back(i);
    }
    selected_n = start;

    // int N = selected_n.size(); 
    int n_max = int(floor(n_total * percentage / 100));

    if (n_max > n_total){n_max = n_total;}

    vec2D newArray;
    for (long unsigned int i = 0; i < selected_n.size(); i++) {
        newArray.push_back(matrix[selected_n[i]]);
    }
    Matrix selection(newArray);

    vector selected_condensed = selection.Sum(COL);

    if (metric == Metric::MSD) {
        Matrix sq_selection = selection.pow(2);
        vector sq_selection_condensed = sq_selection.Sum(AXIS::COLUMN);
        while (static_cast<int>(selected_n.size()) < n_max){
            //select_from_n = np.delete(total_indices, selected_n)
            // new_index_n = get_new_index_n(matrix, metric=metric, selected_condensed,
            //                               sq_selected_condensed, N, 
            //                               select_from_n, n_atoms)
            // sq_selected_condensed += matrix[new_index_n] ** 2
        }
    } else { 
        while (static_cast<int>(selected_n.size()) < n_max){
            // new_index_n = get_new_index_n(matrix, metric, selected_condensed, 
            //                               N, select_from_n)
        }
    }
    // selected_condensed += matrix[new_index_n]
    // selected_n.append(new_index_n)
    // n = len(selected_n)
    
    // return selected_n;

    return total_indices; // Placeholder return
}