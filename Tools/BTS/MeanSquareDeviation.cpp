#include "MeanSquareDeviation.h"

/*
Mean square deviation (MSD) calculation for n-ary objects.
    - Assumes that for matrix NxM that N=M
    - Assumes vector is non-empty

Parameters
----------
matrix : Matrix
    Data matrix.
N_atoms : int
    Number of atoms in the system.

Returns
-------
float
    normalized MSD value.
*/
float MeanSquareDeviation(Matrix matrix, int n_atoms){

    float msd;
    
    // MSD before dividing by N^2
    float sum = 0;

    // Square the argument matrix to get squared matrix
    Matrix sq_matrix = arma::pow(matrix, 2);
    
    // Summate the columns of the matrix and its squared matrix
    rvector c_sum = arma::sum(matrix, COL);
    rvector sq_sum = arma::sum(sq_matrix, COL);

    // N represents number of rows
    int N = matrix.n_rows;

    // Perform summation component of MSD
    for (uword i = 0; i < c_sum.size(); i++) {
        sum += 2 * (N * sq_sum(i) - pow(c_sum(i),2));
    }

    // Calculate non-normalized msd
    msd = sum / pow(N,2);

    // Return normalized MSD value
    return msd / n_atoms;
}


/* Condensed version of Mean square deviation (MSD).

Parameters
----------
c_sum : vector of size n_features
    Column sum of the data. 
sq_sum : vector of size n_features
    Column sum of the squared data.
N : int
    Number of data points.
n_atoms : int
    Number of atoms in the system.

Returns
-------
float
    normalized MSD value.
*/
float MSDCondensed(rvector c_sum, rvector sq_sum, int N, int n_atoms){
    float msd = arma::sum(2 * ((N * sq_sum) - arma::pow(c_sum, 2))) / pow(N,2);
    // float sum = 0;
    // for (uword i = 0; i < c_sum.size(); i++) {
    //     sum += 2 * (N * sq_sum(i) - pow(c_sum(i),2));
    // }
    // float msd = sum / pow(N,2);
    return (msd / n_atoms);
}