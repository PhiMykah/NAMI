#ifndef ESIM_MODULES_H
#define ESIM_MODULES_H
#include <map>
#include <string>
#include "../../../Datatypes/DataContainers.h"

namespace ESIM{
    struct Counters
    {   
        Counters(float a_, float w_a_, float d_, float w_d_, float total_sim_, float total_w_sim_, 
                float total_dis_, float total_w_dis_, float p_, float w_p_)
        {
            a = a_; w_a = w_a_;
            d = d_; w_d = w_d_;
            total_sim = total_sim_; total_w_sim = total_w_sim_;
            total_dis = total_dis_; total_w_dis = total_w_dis_;
            p = p_; w_p = w_p_;
        }

        float a;
        float w_a;
        float d;
        float w_d;
        float total_sim;
        float total_w_sim;
        float total_dis;
        float total_w_dis;
        float p;
        float w_p;
    };

    struct Indices 
    {
        Indices(float bub_nw_, float fai_nw_, float gle_nw_, float ja_nw_, float jt_nw_,
                float rt_nw_, float rr_nw_, float sm_nw_, float ss1_nw_, float ss2_nw_)
        {
            bub_nw = bub_nw_; fai_nw = fai_nw_;
            gle_nw = gle_nw_; ja_nw = ja_nw_;
            jt_nw = jt_nw_; rt_nw = rt_nw_;
            rr_nw = rr_nw_; sm_nw = sm_nw_;
            ss1_nw = ss1_nw_; ss2_nw = ss2_nw_;
        }

        float bub_nw;
        float fai_nw;
        float gle_nw;
        float ja_nw;
        float jt_nw;
        float rt_nw;
        float rr_nw;
        float sm_nw;
        float ss1_nw;
        float ss2_nw;
    };

    // Calculate 1-similarity, 0-similarity, and dissimilarity counters
    Counters CalculateCounters(Matrix c_total, int n_objects, float c_threshold, int w_factor = (int) WFactor::FRACTION);

    // Generate a dict (string->float) map with the similarity indices
    Indices GenSimIndices(Matrix c_total, int n_objects, float c_threshold, int w_factor = (int) WFactor::FRACTION);
}

enum class THRESHOLD {MIN = -2, DISSIMILAR=-1, NONE=0};

// ******************************************
// * Similarity and dissimilarity functions *
// ******************************************

Matrix F_S_Default(Matrix d, int w_factor, int n_objects);

Matrix F_D_Default(Matrix d, int w_factor, int n_objects);

Matrix F_S_Fraction(Matrix d, int w_factor, int n_objects);

Matrix F_D_Fraction(Matrix d, int w_factor, int n_objects);

Matrix F_S_Power(Matrix d, int w_factor, int n_objects);

// Unsure if order of operations with modulo is correct - φ
Matrix F_D_Power(Matrix d, int w_factor, int n_objects);

#endif // !ESIM_MODULES_H