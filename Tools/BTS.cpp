#include "BTS.h"
#include "Esim_Modules.h"

AXIS COL = AXIS::COLUMN; 

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
    // MSD before divinging by N^2
    float sum = 0;

    // Square the argument matrix to get squared matrix
    Matrix sq_matrix = matrix.pow(2);
    
    // Summate the columns of the matrix and its squared matrix
    vector c_sum = matrix.Sum(COL);
    vector sq_sum = sq_matrix.Sum(COL);

    // N represents number of rows
    int N = matrix.N;

    // Perform summation component of MSD
    for (int i = 0; i < c_sum.size(); i++) {
        sum += 2 * (N * sq_sum[i] - pow(c_sum[i],2));
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
float MSDCondensed(vector c_sum, vector sq_sum, int N, int n_atoms){
    float sum = 0;
    for (int i = 0; i < c_sum.size(); i++) {
        sum += 2 * (N * sq_sum[i] - pow(c_sum[i],2));
    }
    float msd = sum / pow(N,2);
    return (msd / n_atoms);
}


/* Calculate the extended comparison of a dataset.
Parameters
----------
matrix : Matrix
    Input data matrix.
metric : {'MSD', 'BUB', 'Fai', 'Gle', 'Ja', 'JT', 'RT', 'RR', 'SM', 'SS1', 'SS2'}
    Metric to use for the extended comparison. Defaults to 'MSD'.
    Available metrics:
    Mean square deviation (MSD), Bhattacharyya's U coefficient (BUB),
    Faiman's coefficient (Fai), Gleason's coefficient (Gle),
    Jaccard's coefficient (Ja), Jaccard-Tanimoto coefficient (JT),
    Rogers-Tanimoto coefficient (RT), Russell-Rao coefficient (RR),
    Simpson's coefficient (SM), Sokal-Sneath 1 coefficient (SS1),
    Sokal-Sneath 2 coefficient (SS2).
N : int, optional
    Number of data points. Defaults to None.
N_atoms : int, optional
    Number of atoms in the system. Defaults to 1.
c_threshold : int, optional
    Coincidence threshold. Defaults to None.
w_factor : {'fraction', 'power_n'}, optional
    Type of weight function that will be used. Defaults to 'fraction'.
    See `esim_modules.calculate_counters` for more information.

Returns 
-------
float
    Extended comparison value.
*/
float ExtendedComparison(
    Matrix matrix, Metric metric, int N, int n_atoms,
    float c_threshold, WFactor w_factor){
    
    // Column sum
    vector c_sum;

    // Calculate the column sum of the matrix
    c_sum = matrix.Sum(COL);

    // Set the number of rows if not provided
    if (N == 0){
        N = matrix.N;
    }

    if (metric == Metric::MSD) {
        // Squared Matrix
        Matrix sq_matrix;
        // Squared matrix column sum
        vector sq_sum;

        sq_matrix = matrix.pow(2);

        sq_sum = sq_matrix.Sum(COL);
        return MSDCondensed(c_sum,sq_sum, N, n_atoms);
    } else {
        Indices esim_dict = GenSimCounters(Matrix(vec2D {c_sum}), N, c_threshold, (int) w_factor);

        switch (metric)
        {
            case Metric::BUB:
                return 1 - esim_dict.bub_nw;
            case Metric::FAI:
                return 1 - esim_dict.fai_nw;
            case Metric::GLE:
                return 1 - esim_dict.gle_nw;
            case Metric::JA:
                return 1 - esim_dict.ja_nw;
            case Metric::JT:
                return 1 - esim_dict.jt_nw;
            case Metric::RT:
                return 1 - esim_dict.rt_nw;
            case Metric::RR:
                return 1 - esim_dict.rr_nw;
            case Metric::SM:
                return 1 - esim_dict.sm_nw;
            case Metric::SS1:
                return 1 - esim_dict.ss1_nw;
            case Metric::SS2:
                return 1 - esim_dict.ss2_nw;
            default:
                break;
        }
        return 1.0;
    }
}


// Calculate the extended comparison of the column sum dataset
float ExtendedComparison(
    vector c_sum, Metric metric, 
    int N, int n_atoms, float c_threshold, 
    WFactor w_factor)
{
    Matrix c_sum_Matrix(vec2D {c_sum});
    
    Indices esim_dict = GenSimCounters(c_sum_Matrix, N, c_threshold, (int) w_factor);

    switch (metric)
    {
        case Metric::BUB:
            return 1 - esim_dict.bub_nw;
        case Metric::FAI:
            return 1 - esim_dict.fai_nw;
        case Metric::GLE:
            return 1 - esim_dict.gle_nw;
        case Metric::JA:
            return 1 - esim_dict.ja_nw;
        case Metric::JT:
            return 1 - esim_dict.jt_nw;
        case Metric::RT:
            return 1 - esim_dict.rt_nw;
        case Metric::RR:
            return 1 - esim_dict.rr_nw;
        case Metric::SM:
            return 1 - esim_dict.sm_nw;
        case Metric::SS1:
            return 1 - esim_dict.ss1_nw;
        case Metric::SS2:
            return 1 - esim_dict.ss2_nw;
        default:
            break;
    }
    return 1.0;
}

// Calculate the extended comparison of the column sum and square column sum of datasets
float ExtendedComparison(
    vector c_sum, vector sq_sum, 
    Metric metric, int N, int n_atoms,
    float c_threshold, WFactor w_factor)
{
    // Only use sq_sum if the metric is MSD
    if (metric == Metric::MSD) {
        return MSDCondensed(c_sum, sq_sum, N, n_atoms);
    }    

    return ExtendedComparison(c_sum, metric, N, n_atoms, c_threshold, w_factor);
}

// Complementary similarity is calculating the similarity of a set 
// without one object or observation using metrics in the extended comparison.
// The greater the complementary similarity, the more representative the object is.
Matrix CalculateCompSim(Matrix matrix, Metric metric, int n_atoms){
    if (metric == Metric::MSD) {
        return CSimMSD(matrix, n_atoms);
    }

    int N = matrix.N;

    Matrix sq_data_total = matrix.pow(2);

    vector c_sum_total = matrix.Sum(AXIS::COLUMN);

    vector sq_sum_total = sq_data_total.Sum(AXIS::COLUMN); 

    vector values;
    for (int row = 0; row < N; row++){
        float value = ExtendedComparison(c_sum_total - matrix[row], metric, N-1, n_atoms);
        values.push_back(value);
    }

    return Matrix(vec2D {values});
}

// Simplified complementary similarity calculation if metric is MSD
Matrix CSimMSD(Matrix matrix, int n_atoms){
    int N = matrix.N;
    Matrix sq_data = matrix.pow(2);

    vector c_sum = matrix.Sum(AXIS::COLUMN);
    vector sq_sum = sq_data.Sum(AXIS::COLUMN);

    Matrix comp_csum = c_sum - matrix;
    Matrix comp_sqsum = sq_sum - sq_data;

    Matrix total = vec2D{(2 * ((N-1) * comp_sqsum - comp_csum.pow(2))).Sum(AXIS::ROW)};

    Matrix comp_msd = total / std::pow(N-1, 2);

    Matrix norm_msd = comp_msd / n_atoms;

    return norm_msd;
}


// Calculates the medoid of a dataset using the metrics in extended comparison.
// Medoid is the most representative object of a set.
int CalculateMedoid(Matrix matrix, Metric metric, int n_atoms){
    if (metric == Metric::MSD) {
        // Returns the indicies where the maximum value occurs
        vector csim = vector(CSimMSD(matrix, n_atoms));
        vector::iterator it = max_element(csim.begin(), csim.end());
        int index = std::distance(csim.begin(), it);
        return index;
    }
    int N = matrix.N;
    Matrix sq_data_total = matrix.pow(2);
    vector c_sum_total = matrix.Sum(AXIS::COLUMN);
    vector sq_sum_total = sq_data_total.Sum(AXIS::COLUMN);
    int index = N + 1;
    int max_dissim = -1;

    for (int row = 0; row < N; row++){
        float value = ExtendedComparison(c_sum_total - matrix[row], metric, N-1, n_atoms);
        if (value > max_dissim) {
            max_dissim = value;
            index = row;
        }
    }

    return index;
}

// Calculates the outliers of a dataset using the metrics in extended comparison.
// Outliers are the least representative objects of a set.
int CalculateOutlier(Matrix matrix, Metric metric, int n_atoms){
    if (metric == Metric::MSD) {
        // Returns the indicies where the minimum value occurs
        vector csim = vector(CSimMSD(matrix, n_atoms));
        vector::iterator it = min_element(csim.begin(), csim.end());
        int index = std::distance(csim.begin(), it);
        return index;
    }
    int N = matrix.N;
    Matrix sq_data_total = matrix.pow(2);
    vector c_sum_total = matrix.Sum(AXIS::COLUMN);
    vector sq_sum_total = sq_data_total.Sum(AXIS::COLUMN);
    int index = N + 1;
    int max_dissim = INT_MAX;

    for (int row = 0; row < N; row++){
        float value = ExtendedComparison(c_sum_total - matrix[row], metric, N-1, n_atoms);
        if (value < max_dissim) {
            max_dissim = value;
            index = row;
        }
    }

    return index;

}

/*
// Trims a desired percentage of outliers (most dissimilar) from the dataset 
// by calculating largest complement similarity.
Matrix TrimOutliers(
    Matrix matrix, float n_trimmed, Metric metric, 
    int n_atoms=1, Criterion criterion = Criterion::COMP_SIM){

}

// Overloaded version of default TrimOutliers
Matrix TrimOutliers(
    Matrix matrix, int n_trimmed, Metric metric,
    int n_atoms=1, Criterion criterion = Criterion::COMP_SIM){

}

// Selects a diverse subset of the data using the complementary similarity.
std::vector<int> DiversitySelection(
    Matrix matrix, int percentage, Metric metric,
    DiversitySeed start = DiversitySeed::MEDOID, int n_atoms = 1){

}

// Function to get the new index to add to the selected indices
int GetNewIndexN(Matrix matrix, Metric metric, Matrix select_condensed,
    int N, std::vector<int> select_from_n, int n_atoms = 1){

}

// Overloaded version of default GetNewIndexN
int GetNewIndexN(Matrix matrix, Metric metric, Matrix select_condensed,
    int N, std::vector<int> select_from_n, Matrix sq_selected_condensed,
    int N_atoms = 1){

}

// Aligns trajectory using uniform or kronecker alignment.
Matrix AlignTraj(Matrix data, int n_atoms, AlignMethod align_method = AlignMethod::UNI){

}

// Aligns the frames in the trajectory to the reference frame.
Matrix EquilAlign(
    std::vector<int> indices, int sieve, 
    std::string input_top, std::string input_traj,
    std::string mdana_atomsel, std::string ccptraj_atomsel,
    int ref_index){

}
*/