#ifndef NEW_INDEX_H
#define NEW_INDEX_H
#include "BTS.h"
#include "ExtendedComparison.h"

// Function to get the new index to add to the selected indices
uword GetNewIndexN(Matrix matrix, Metric metric, rvector select_condensed,
    uword N, index_vec select_from_n, int n_atoms = 1);
uword GetNewIndexN(Matrix matrix, Metric metric, rvector select_condensed,
    rvector sq_selected_condensed, uword N, index_vec select_from_n,
    int n_atoms = 1);

#endif // !NEW_INDEX_H