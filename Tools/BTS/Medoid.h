#ifndef MEDOID_H
#define MEDOID_H
#include "BTS.h"
#include "ExtendedComparison.h"
#include "ComplementarySimilarity.h"
#include <stdio.h>

// Calculates the medoid of a dataset using the metrics in extended comparison.
// Medoid is the most representative object of a set.
int CalculateMedoid(Matrix matrix, Metric metric, int n_atoms = 1);

#endif // !MEDOID_H