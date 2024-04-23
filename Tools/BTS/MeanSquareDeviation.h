#ifndef MEAN_SQUARE_DEVIATION_H
#define MEAN_SQUARE_DEVIATION_H
#include "BTS.h"

// Mean square deviation (MSD) calculation for n-ary objects.
float MeanSquareDeviation(Matrix matrix, int n_atoms);

// Condensed version of Mean square deviation (MSD).
float MSDCondensed(rvector c_sum, rvector sq_sum, int N, int n_atoms);
#endif // !MEAN_SQUARE_DEVIATION_H