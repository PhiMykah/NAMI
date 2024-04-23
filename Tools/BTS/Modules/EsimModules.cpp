#include "EsimModules.h"

/*
Calculate 1-similarity, 0-similarity, and dissimilarity counters

Parameters
---------
c_total : Matrix
    Matrix (n_objects, n_features) containing the sums of each column of the fingerprint matrix.

n_objects : int
    Number of objects to be compared.

c_threshold : float
    Coincidence threshold.
    Treated as percentage if in range (0,1)
    Treated as integer if greater than 1
    THRESHOLD::None or THRESHOLD::MIN : Default, c_threshold = n_objects % 2
    THRESHOLD::DISSIMILAR : c_threshold = ceil(n_objects / 2)
    Otherwise : Integer number < n_objects

w_factor : int ;
    Type of weight function that will be used.
    W_FACTOR::FRACTION = -1 : similarity = d[k]/n
                    dissimilarity = 1 - (d[k] - n_objects % 2)/n_objects
    int : similarity = n**-(n_objects - d[k])
                dissimilarity = n**-(d[k] - n_objects % 2)
    Otherwise : similarity = dissimilarity = 1

Returns
-------
counters : Counters
    Struct object with the weighted and non-weighted counters.
*/
Counters CalculateCounters(Matrix c_total, int n_objects, float c_threshold, int w_factor){
    if ((0 < c_threshold) && (c_threshold < 1)) { c_threshold *= n_objects; } 
    else {
        switch ((int) c_threshold)
        {
        // Set threshold to even remainder if NONE or MIN is selected
        case (int) THRESHOLD::NONE:
        case (int) THRESHOLD::MIN:
            c_threshold = n_objects % 2;
            break;
        // Set threshold to half the number of objects, rounded up
        case (int) THRESHOLD::DISSIMILAR:
            c_threshold = std::ceil(n_objects / 2);
            break;
        default:
            // Ensure that the threshold is less than n_objects
            if (c_threshold >= n_objects) {
                c_threshold = n_objects-1;
            }
            break;
        }
    }

    // Implementation of f_s and f_d emulates python mdance, may be 
    // overall inefficient for c++ - φ
    Matrix (*f_s)(Matrix, int, int); 
    Matrix (*f_d)(Matrix, int, int);
    
    switch (w_factor)
    {
    case 0:
        f_s = &F_S_Default;
        f_d = &F_D_Default;
        break;
    case -1:
        f_s = &F_S_Fraction;
        f_d = &F_D_Fraction;
        break;
    default:
        f_s = &F_S_Power;
        f_d = &F_D_Power;
        break;
    }
    
    /*
    This section implements MDANCE's logic almost identically, but with a different datatype
    */

    index_vec a_indices = arma::find((2 * c_total - n_objects) > c_threshold);
    index_vec d_indices = arma::find((n_objects - 2 * c_total) > c_threshold);
    index_vec dis_indices = arma::find((arma::abs(2 * c_total - n_objects) <= c_threshold));

    // Calculate overall sum by summing all of the columns and summing all of the remaining rows
    float a = a_indices.n_elem;
    float d = d_indices.n_elem;
    float total_dis = dis_indices.n_elem;
    
    Matrix a_w_array = f_s(2 * c_total.elem(a_indices) - n_objects, w_factor, n_objects);
    Matrix d_w_array = f_s(arma::abs(2 * c_total.elem(d_indices) - n_objects), w_factor, n_objects);
    Matrix total_w_dis_array = f_d(arma::abs(2 * c_total.elem(dis_indices) - n_objects), w_factor, n_objects);

    float w_a = MatSum(a_w_array);
    float w_d = MatSum(d_w_array);
    float total_w_dis = MatSum(total_w_dis_array);

    float total_sim = a + d;
    float total_w_sim = w_a + w_d;
    float p = total_sim + total_dis;
    float w_p = total_w_sim + total_w_dis;

    // Return counters struct
    Counters counters(
        a, w_a, 
        d, w_d, 
        total_sim, total_w_sim,
        total_dis, total_w_dis, 
        p, w_p
    );
    return counters;
}

/*
Generate a struct (attribute -> float) with the similarity indices

Parameters
----------
c_total : Matrix
    Matrix (n_objects, n_features) containing the sums of each column of the fingerprint matrix.
n_objects : int
    Number of objects to be compared.
c_threshold : float/int
    Coincidence threshold.
w_factor : int
    Type of weight function that will be used.

Returns
-------
counters : Counters
    Counters struct with the similarity indices.

Notes
-----
Available indices:
BUB: Baroni-Urbani-Buser, Fai: Faith, Gle: Gleason, Ja: Jaccard,
JT: Jaccard-Tanimoto, RT: Rogers-Tanimoto, RR: Russel-Rao
SM: Sokal-Michener, SSn: Sokal-Sneath n
*/
Indices GenSimCounters(Matrix c_total, int n_objects, float c_threshold, int w_factor){
    Counters counters = CalculateCounters(c_total, n_objects, c_threshold, w_factor);

    float bub_nw = (std::pow((counters.w_a * counters.w_d), 0.5) + counters.w_a)/
             (std::pow((counters.a * counters.d), 0.5) + counters.a + counters.total_dis);
    float fai_nw = (counters.w_a + 0.5 * counters.w_d)/
             (counters.p);
    float gle_nw = (2 * counters.w_a)/
             (2 * counters.a + counters.total_dis);
    float ja_nw = (3 * counters.w_a)/
            (3 * counters.a + counters.total_dis);
    float jt_nw = (counters.w_a)/
            (counters.a + counters.total_dis);
    float rt_nw = (counters.total_w_sim)/
            (counters.p + counters.total_dis);
    float rr_nw = (counters.w_a)/
            (counters.p);
    float sm_nw = (counters.total_w_sim)/
            (counters.p);
    float ss1_nw = (counters.w_a)/
             (counters.a + 2 * counters.total_dis);
    float ss2_nw = (2 * counters.total_w_sim)/
             (counters.p + counters.total_sim);

    Indices indices(
        bub_nw, fai_nw, gle_nw,
        ja_nw, jt_nw, rt_nw,
        rr_nw, sm_nw, ss1_nw, ss2_nw
    );
    return indices;
}

// ******************************************
// * Similarity and dissimilarity functions *
// ******************************************

Matrix F_S_Default(Matrix d, int w_factor, int n_objects){return Matrix(1,1, arma::fill::ones);}

Matrix F_D_Default(Matrix d, int w_factor, int n_objects){return Matrix(1,1, arma::fill::ones);}

Matrix F_S_Fraction(Matrix d, int w_factor, int n_objects){return d/n_objects;}

// Unsure if order of operations with modulo is correct - φ
Matrix F_D_Fraction(Matrix d, int w_factor, int n_objects){return 1 - (d - n_objects % 2)/n_objects;}

Matrix F_S_Power(Matrix d, int w_factor, int n_objects){
    Matrix w_fac(d.n_rows, d.n_cols, arma::fill::value(w_factor)); 
    return arma::pow(w_fac,-1 * (n_objects-d));}

// Unsure if order of operations with modulo is correct - φ
Matrix F_D_Power(Matrix d, int w_factor, int n_objects){
    Matrix w_fac(d.n_rows, d.n_cols, arma::fill::value(w_factor));
    return arma::pow(w_fac,-1 * (d - (n_objects % 2)));}