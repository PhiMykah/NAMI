/*
                            iSIM_MODULES
    ----------------------------------------------------------------------
    
    Miranda-Quintana Group, Department of Chemistry, University of Florida 
    
    ----------------------------------------------------------------------
    
    Please, cite the original paper on iSIM:
*/
#ifndef ISIM_MODULES_H
#define ISIM_MODULES_H
#include <map>
#include <string>
#include "../../../Datatypes/DataContainers.h"

namespace ISIM{
    struct Counters
    {   
        Counters(float a_, float d_, float total_sim_, 
                float total_dis_, float p_)
        {
            a = a_; d = d_;
            total_sim = total_sim_;
            total_dis = total_dis_;
            p = p_;
        }

        float a;
        float d;
        float total_sim;
        float total_dis;
        float p;
    };

    struct Indices 
    {
        Indices(float ac_, float bub_, float fai_, float gle_, float ja_,
                float jt_, float rt_, float rr_, float sm_, float ss1_, float ss2_)
        {
            ac = ac_; bub = bub_;
            fai = fai_; gle = gle_;
            ja = ja_; jt = jt_;
            rt = rt_; rr = rr_;
            sm = sm_; ss1 = ss1_;
            ss2 = ss2_;
        }

        float ac;
        float bub;
        float fai;
        float gle;
        float ja;
        float jt;
        float rt;
        float rr;
        float sm;
        float ss1;
        float ss2;
    };

    Counters CalculateCounters(Matrix data, int k = 1);

    Counters CalculateCounters(vector c_total, int n_objects, int k = 1);

    Indices GenSimIndices(Matrix data, int n_objects = 0, int k = 1);
}

float CalculateIsim(Matrix data, Metric n_ary = Metric::RR);

float CalculateIsim(vector c_total, int n_objects, Metric n_ary = Metric::RR);

#endif // !ISIM_MODULES_H