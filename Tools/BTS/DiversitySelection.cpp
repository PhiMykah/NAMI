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
index_vec DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    DiversitySeed start, int n_atoms)
{
    index_vec selected_n(1);
    uword n_total = matrix.n_rows;
    int seed;

    switch (start)
    {
    case DiversitySeed::OUTLIER:
        seed = CalculateOutlier(matrix, metric, n_atoms);
        break;
    case DiversitySeed::RANDOM:
        seed = rand() % n_total;
    default:
        seed = CalculateMedoid(matrix, metric, n_atoms);
        break;
    }

    if ((seed > 0) && ((uword) seed >= n_total)) {
        fprintf(stderr, "selected_n value %i for metric %s is out of range [%i, %llu), defaulting to %llu\n", 
                seed, toStr(metric).c_str(), 0, n_total, n_total-1);
        seed = n_total - 1;
    }
    if (seed < 0) {
        fprintf(stderr, "selected_n value %i for metric %s is out of range [%i, %llu), defaulting to %i\n", 
                seed, toStr(metric).c_str(), 0, n_total, 0);
        seed = 0;
    }

    selected_n(0) = seed;

    return DiversitySelection(matrix, percentage, metric, selected_n, n_atoms);
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
start : std::vector<int>
    Seed vector of diversity selection.
N_atoms : int, optional
    Number of atoms in the system. Defaults to 1.

Returns
-------
list
    List of indices of the selected data.
*/
index_vec DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    index_vec start, int n_atoms)
{   
    // Variable declarations 
    index_vec selected_n = start;
    uword new_index_n;
    rvector sq_selection_condensed;
    uword n_total = matrix.n_rows;
    uword prev_size;
    index_vec total_indices = arma::regspace<index_vec>(0, n_total-1);

    uword N = selected_n.size(); 
    uword n_max = (uword)floor(n_total * percentage / 100);

    if (n_max > n_total){n_max = n_total;}

    Matrix selection(selected_n.size(), matrix.n_cols, arma::fill::zeros);

    for (long unsigned int i = 0; i < selected_n.size(); i++) {
        selection.row(i) = matrix.row(selected_n[i]);
    }

    rvector selected_condensed = arma::sum(selection,COL);

    if (metric == Metric::MSD) {
            sq_selection_condensed = arma::sum(arma::pow(selection,2),COL);
    }
    while (selected_n.size() < n_max){
        //select_from_n = np.delete(total_indices, selected_n)
        //Removing rows each time leads to issues 
        index_vec select_from_n(total_indices);
        select_from_n.shed_rows(selected_n); 

        if (metric == Metric::MSD){
            // new_index_n = get_new_index_n(matrix, metric=metric, selected_condensed,
            //                               sq_selected_condensed, N, 
            //                               select_from_n, n_atoms)
            new_index_n = GetNewIndexN(
                matrix, metric, selected_condensed, 
                sq_selection_condensed, N, select_from_n, n_atoms);

            // sq_selected_condensed += matrix[new_index_n] ** 2
            sq_selection_condensed = sq_selection_condensed + pow(matrix.row(new_index_n), 2);
        } else {
            // new_index_n = get_new_index_n(matrix, metric, selected_condensed, 
            //                               N, select_from_n)
            new_index_n = GetNewIndexN(
                matrix, metric, selected_condensed, 
                N, select_from_n, n_atoms);

        }
        // selected_condensed += matrix[new_index_n]
        selected_condensed = selected_condensed + matrix.row(new_index_n);

        // selected_n.append(new_index_n)
        prev_size = selected_n.size();
        selected_n.resize(prev_size+1);
        selected_n(prev_size) = new_index_n;

        // n = len(selected_n)
        N = selected_n.size();
    }

    
    return selected_n;
}