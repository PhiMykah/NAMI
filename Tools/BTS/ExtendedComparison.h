#ifndef EXTENDED_COMPARISON_H
#define EXTENDED_COMPARISON_H
#include "BTS.h"
#include "MeanSquareDeviation.h"

// Calculate the extended comparison of a dataset. 
float ExtendedComparison(
    Matrix matrix, Metric metric = Metric::MSD, int N = 0, int n_atoms = 1,
    float c_threshold = 0, WFactor w_factor = WFactor::FRACTION);

// Calculate the extended comparison of the column sum dataset
float ExtendedComparison(
    vector c_sum, Metric metric = Metric::MSD, 
    int N = 0, int n_atoms = 1, float c_threshold = 0, 
    WFactor w_factor = WFactor::FRACTION);

// Calculate the extended comparison of the column sum and square column sum of datasets
float ExtendedComparison(
    vector c_sum, vector sq_sum, 
    Metric metric = Metric::MSD, int N = 0, int n_atoms = 1,
    float c_threshold = 0, WFactor w_factor = WFactor::FRACTION);

#endif // !EXTENDED_COMPARISON_H