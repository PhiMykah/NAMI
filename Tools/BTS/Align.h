#ifndef ALIGN_H
#define ALIGN_H
#include "BTS.h"

// Aligns trajectory using uniform or kronecker alignment.
Matrix AlignTraj(Matrix data, int n_atoms, AlignMethod align_method = AlignMethod::UNI);

// Aligns the frames in the trajectory to the reference frame.
Matrix EquilAlign(
    std::vector<int> indices, int sieve, 
    std::string input_top, std::string input_traj,
    std::string mdana_atomsel, std::string ccptraj_atomsel,
    int ref_index);
    
#endif // !ALIGN_H