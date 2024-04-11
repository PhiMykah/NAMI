#include "ExtendedComparison.h"

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
c_threshold : float, optional
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


/* Calculate the extended comparison of the column sum dataset

Parameters
----------
c_sum : vector
    Input column sum vector.
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
c_threshold : float, optional
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


/* Calculate the extended comparison of the column sum and square column sum of datasets

Parameters
----------
c_sum : vector
    Input column sum vector.
sq_sum : vector
    Input square column sum vector.
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
c_threshold : float, optional
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