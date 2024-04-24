#include "IsimModules.h"

/*
Calculate 1-similarity, 0-similarity, and dissimilarity counters

Parameters
----------
data : Matrix
    Matrix where dach vector contains the binary object 

k : int
    Integer indicating the 1/k power used to approximate the average of the
    similarity values elevated to 1/k.

Returns
-------
counters : ISIM::Counters
    Struct object with the weighted and non-weighted counters.

*/
ISIM::Counters ISIM::CalculateCounters(Matrix data, int k)
{
    vector c_total = arma::sum(data, (int)AXIS::COLUMN);
    int n_objects = data.n_rows;

    return ISIM::CalculateCounters(c_total, n_objects, k);
}

/*
Calculate 1-similarity, 0-similarity, and dissimilarity counters

Parameters
----------
c_toal : vector
    Vector with the columnwise sum, if so specify n_objects

n_objects : int
    Number of objects.

k : int
    Integer indicating the 1/k power used to approximate the average of the
    similarity values elevated to 1/k.

Returns
-------
counters : ISIM::Counters
    Struct object with the weighted and non-weighted counters.

*/
ISIM::Counters ISIM::CalculateCounters(vector c_total, int n_objects, int k)
{
    vector a_array = c_total * (c_total - 1) / 2;
    vector off_coincidences = n_objects - c_total;
    vector d_array = off_coincidences * (off_coincidences -1) / 2;
    vector dis_array = off_coincidences * c_total;

    float a = MatSum(arma::pow(a_array, 1/k));
    float d = MatSum(arma::pow(d_array, 1/k));
    float total_dis = MatSum(arma::pow(dis_array, 1/k));

    float total_sim = a + d;
    float p = total_sim + total_dis;
    
    return ISIM::Counters(a, d, total_sim, total_dis, p);
}

/* 
Calculate the iSIM index for RR, JT, or SM

Parameters
----------
data : Matrix
    Each vector contains the binary object.

n_ary : Metric
    Desired similarity index to calculate the iSIM from. 
    Only RR, JT, or SM are available. For other indexes use gen_sim_dict.

Returns
-------
isim : float
    iSIM index for the specified similarity index.
*/
float CalculateIsim(Matrix data, Metric n_ary)
{
    vector c_total = arma::sum(data, (int)AXIS::COLUMN);
    int n_objects = data.n_rows;
    return CalculateIsim(c_total, n_objects, n_ary);
}

/*
Calculate the iSIM index for RR, JT, or SM

Parameters
----------
data : vector
    vector of the columnwise sum

n_objects : int
    Number of objects.

n_ary : Metric
    Desired similarity index to calculate the iSIM from. 
    Only RR, JT, or SM are available. For other indexes use gen_sim_dict.

Returns
-------
isim : float
    iSIM index for the specified similarity index.
*/

float CalculateIsim(vector c_total, int n_objects, Metric n_ary)
{
    float a = MatSum(c_total * (c_total - 1) / 2);
    float p = n_objects * (n_objects - 1) * c_total.size() / 2;
    if (n_ary == Metric::JT)
    {
        vector off_coincidences = n_objects - c_total;
        float total_dis = MatSum(off_coincidences * c_total);

        return a/(a + total_dis);

    } else if (n_ary == Metric::SM)
    {
        vector off_coincidences = n_objects - c_total;
        float d = MatSum(off_coincidences * (off_coincidences - 1) / 2);
        
        return (a + d)/p;
    }
    if (n_ary != Metric::RR) {
        fprintf(stderr, "Unexpected Metric \'%s\', defaulting to \'%s\'",
        toStr(n_ary).c_str(), toStr(Metric::RR).c_str());
    } 

    return a/p;
}


/*
Calculate a dictionary containing all the available similarity indexes

Parameters
----------
See CalculateCounters.

Returns
-------
sim_dict : ISIM::Indices
    Struct with the weighted and non-weighted similarity indexes.
*/
ISIM::Indices ISIM::GenSimIndices(Matrix data, int n_objects, int k)
{  
    ISIM::Counters counters = ISIM::CalculateCounters(data, n_objects, k);
    
    float ac = (2/M_PI) * std::asin(std::sqrt(counters.total_sim/counters.p));
    float bub = ((std::sqrt(counters.a * counters.d) + counters.a)/\
                 (std::sqrt(counters.a * counters.d)) + counters.a + counters.total_dis);
    float fai = (counters.a + 0.5 * counters.d)/\
                (counters.p);
    float gle = (2 * counters.a)/\
          (2 * counters.a + counters.total_dis);
    float ja = (3 * counters.a)/\
         (3 * counters.a + counters.total_dis);
    float jt = (counters.a)/\
         (counters.a + counters.total_dis);
    float rt = (counters.total_sim)/\
         (counters.p + counters.total_dis);
    float rr = (counters.a)/\
         (counters.p);
    float sm = (counters.total_sim)/\
         (counters.p);
    float ss1 = (counters.a)/\
          (counters.a + 2 * counters.total_dis);
    float ss2 = (2 * counters.total_sim)/\
          (counters.p + counters.total_sim);


    // Dictionary with all the results
    return ISIM::Indices(ac, bub, fai, gle, ja,
                                   jt, rt, rr, sm, ss1, ss2);
    // Indices = {'Fai':fai, 'Gle':gle, 'Ja':ja,
    //            'JT':jt, 'RT':rt, 'RR':rr, 'SM':sm, 'SS1':ss1, 'SS2':ss2}
}
