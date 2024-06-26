#ifndef DIVERSITY_SELECTION_H
#define DIVERSITY_SELECTION_H
#include "BTS.h"
#include "Outlier.h"
#include "Medoid.h"
#include "NewIndex.h"

// Selects a diverse subset of the data using the complementary similarity.
index_vec DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    DiversitySeed start = DiversitySeed::MEDOID, int n_atoms = 1);

// Selects a diverse subset of the data using the complementary similarity.
index_vec DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    index_vec start, int n_atoms = 1);
    
#endif // !DIVERSITY_SELECTION_H