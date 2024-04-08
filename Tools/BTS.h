// BTS: Behind the Scenes tools functions
// File/function names and types would need to be adopted to the cpptraj style
#ifndef BTS_H
#define BTS_H
#include "Data_Containers.h"
#include <stdio.h>
#include <iostream>

// Mean square deviation (MSD) calculation for n-ary objects.
float MeanSquareDeviation(Matrix matrix, int n_atoms);

// Condensed version of Mean square deviation (MSD).
float MSDCondensed(std::vector<float> c_sum, std::vector<float> sq_sum, int N, int n_atoms);

// Calculate the extended comparison of a dataset. 
float ExtendedComparison(
    Matrix matrix, Metric metric = Metric::MSD, int N = 0, int n_atoms = 1,
    float c_threshold = 0, WFactor w_factor = WFactor::FRACTION);

// Calculate the extended comparison of the column sum dataset
float ExtendedComparison(
    std::vector<float> c_sum, Metric metric = Metric::MSD, 
    int N = 0, int n_atoms = 1, float c_threshold = 0, 
    WFactor w_factor = WFactor::FRACTION);

// Calculate the extended comparison of the column sum and square column sum of datasets
float ExtendedComparison(
    std::vector<float> c_sum, std::vector<float> sq_sum, 
    Metric metric = Metric::MSD, int N = 0, int n_atoms = 1,
    float c_threshold = 0, WFactor w_factor = WFactor::FRACTION);

// Complementary similarity is calculating the similarity of a set 
// without one object or observation using metrics in the extended comparison.
// The greater the complementary similarity, the more representative the object is.
Matrix CalculateCompSim(Matrix matrix, Metric metric, int n_atoms = 1);

// Simplified complementary similarity calculation if metric is MSD
Matrix CSimMSD(Matrix matrix, int n_atoms = 1);

// Calculates the medoid of a dataset using the metrics in extended comparison.
// Medoid is the most representative object of a set.
int CalculateMedoid(Matrix matrix, Metric metric, int n_atoms = 1);

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

// Selects a diverse subset of the data using the complementary similarity.
std::vector<int> DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    DiversitySeed start = DiversitySeed::MEDOID, int n_atoms = 1);

// Function to get the new index to add to the selected indices
int GetNewIndexN(Matrix matrix, Metric metric, Matrix select_condensed,
    int N, std::vector<int> select_from_n, int n_atoms = 1);
int GetNewIndexN(Matrix matrix, Metric metric, Matrix select_condensed,
    int N, std::vector<int> select_from_n, Matrix sq_selected_condensed,
    int N_atoms = 1);

// Aligns trajectory using uniform or kronecker alignment.
Matrix AlignTraj(Matrix data, int n_atoms, AlignMethod align_method = AlignMethod::UNI);

// Aligns the frames in the trajectory to the reference frame.
Matrix EquilAlign(
    std::vector<int> indices, int sieve, 
    std::string input_top, std::string input_traj,
    std::string mdana_atomsel, std::string ccptraj_atomsel,
    int ref_index);

#endif // !BTS_H