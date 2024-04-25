#include "DataContainers.h"

std::string toStr(Metric metric)
{   
    switch (metric){

    case Metric::MSD:
        return("MSD");
    case Metric::BUB:
        return("BUB");
    case Metric::FAI:
        return("FAI");
    case Metric::GLE:
        return("GLE");
    case Metric::JA:
        return("JA");
    case Metric::JT:
        return("JT");
    case Metric::RT:
        return("RT");
    case Metric::RR:
        return("RR");
    case Metric::SM:
        return("SM");
    case Metric::SS1:
        return("SS1");
    case Metric::SS2:
        return("SS2");
    default:
        return("Unknown Metric");
    }

    return("Unknown Metric");
}

void sortRows(Matrix mat, float (*key)(rvector v), uword l_index, uword r_index, bool reverse)
{   
    auto greater_than = [](float val, float pivot)->bool{return val > pivot;};
    auto less_than = [](float val, float pivot)->bool{return val < pivot;};
    auto i_comparison = less_than;
    auto j_comparison = greater_than;
    uword i, j, m_index;
    float pivot;
    i = l_index;
    j = r_index;
    m_index = l_index + (r_index - l_index) / 2;
    pivot = key(mat.row(m_index));

    if (reverse == true) {
        auto i_comparison = greater_than;
        auto j_comparison = less_than;
    }

    while(i < r_index || j > l_index){
        while (i_comparison(key(mat.row(i)),pivot)) {
            i++;
        }
        while (j_comparison(key(mat.row(j)),pivot)) {
            j--;
        }

        if (i <= j) {
            mat.swap_rows(i, j);
            i++;
            j--;
        } else {
            if (i < r_index) {
                sortRows(mat, key, i, r_index, reverse);
            }
            if (j > l_index) {
                sortRows(mat, key, l_index, j, reverse);
            }

            return;
        }
    }
}
