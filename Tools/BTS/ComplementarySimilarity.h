#ifndef COMPLEMENTARY_SIMILARITY_H
#define COMPLEMENTARY_SIMILARITY_H
#include "BTS.h"
#include "ExtendedComparison.h"

// Complementary similarity is calculating the similarity of a set 
// without one object or observation using metrics in the extended comparison.
// The greater the complementary similarity, the more representative the object is.
Matrix CalculateCompSim(Matrix matrix, Metric metric, int n_atoms = 1);

// Simplified complementary similarity calculation if metric is MSD
Matrix CSimMSD(Matrix matrix, int n_atoms = 1);

#endif // !COMPLEMENTARY_SIMILARITY_H