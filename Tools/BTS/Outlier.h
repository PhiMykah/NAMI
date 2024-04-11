#ifndef OUTLIER_H
#define OUTLIER_H
#include "BTS.h"
#include "ExtendedComparison.h"
#include "ComplementarySimilarity.h"
#include "Medoid.h"

// Calculates the outliers of a dataset using the metrics in extended comparison.
// Outliers are the least representative objects of a set.
int CalculateOutlier(Matrix matrix, Metric metric, int n_atoms = 1);

// Trims a desired percentage of outliers (most dissimilar) from the dataset 
// by calculating largest complement similarity.
Matrix TrimOutliers(
    Matrix matrix, float percent_trimmed, Metric metric, 
    int n_atoms=1, Criterion criterion = Criterion::COMP_SIM);
Matrix TrimOutliers(
    Matrix matrix, int n_trimmed, Metric metric,
    int n_atoms=1, Criterion criterion = Criterion::COMP_SIM);

#endif // !OUTLIER_H