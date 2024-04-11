#ifndef NEW_INDEX_H
#define NEW_INDEX_H
#include "BTS.h"

// Function to get the new index to add to the selected indices
int GetNewIndexN(Matrix matrix, Metric metric, Matrix select_condensed,
    int N, std::vector<int> select_from_n, int n_atoms = 1);
int GetNewIndexN(Matrix matrix, Metric metric, Matrix select_condensed,
    int N, std::vector<int> select_from_n, Matrix sq_selected_condensed,
    int N_atoms = 1);

#endif // !NEW_INDEX_H