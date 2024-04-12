#include "Align.h"


/* 
Aligns trajectory using uniform or kronecker alignment.

Parameters
----------
data : Matrix
    Input data matrix.
N_atoms : int
    number of atoms in the system
align_method : {'uni', 'kron', 'uniform', 'kronecker'}, optional
    alignment method, by default "uni"

Returns
-------
array-like of shape (n_samples, n_features)
    matrix of aligned data
*/
// Matrix AlignTraj(Matrix data, int n_atoms, AlignMethod align_method){
// }

/*
// Aligns the frames in the trajectory to the reference frame.
Matrix EquilAlign(
    std::vector<int> indices, int sieve, 
    std::string input_top, std::string input_traj,
    std::string mdana_atomsel, std::string ccptraj_atomsel,
    int ref_index){

}
*/