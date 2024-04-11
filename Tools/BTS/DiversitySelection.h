#ifndef DIVERSITY_SELECTION_H
#define DIVERSITY_SELECTION_H
#include "BTS.h"
#include "Outlier.h"
#include "Medoid.h"
#include "NewIndex.h"

// Selects a diverse subset of the data using the complementary similarity.
std::vector<int> DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    DiversitySeed start = DiversitySeed::MEDOID, int n_atoms = 1);

// Selects a diverse subset of the data using the complementary similarity.
std::vector<int> DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    std::vector<int> start, int n_atoms = 1);
    
#endif // !DIVERSITY_SELECTION_H