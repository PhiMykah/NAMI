#ifndef MEAN_SQUARE_DEVIATION_H
#define MEAN_SQUARE_DEVIATION_H
#include "BTS.h"

// Mean square deviation (MSD) calculation for n-ary objects.
float MeanSquareDeviation(Matrix matrix, int n_atoms);

// Condensed version of Mean square deviation (MSD).
float MSDCondensed(std::vector<float> c_sum, std::vector<float> sq_sum, int N, int n_atoms);
#endif // !MEAN_SQUARE_DEVIATION_H