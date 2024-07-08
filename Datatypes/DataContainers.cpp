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
    bool (*i_comparison)(float, float) = less_than;
    bool (*j_comparison)(float, float) = greater_than;
    uword i, j, m_index;
    float pivot;
    i = l_index;
    j = r_index;
    m_index = l_index + (r_index - l_index) / 2;
    pivot = key(mat.row(m_index));

    if (reverse == true) {
        i_comparison = greater_than;
        j_comparison = less_than;
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

double Euclidian(rvector A, rvector B){
    // Distance = sqrt(dot(A, A) - 2 * dot(A, B) + dot(B, B))
    return sqrt(arma::sum(arma::pow(A-B, 2)));
}

double Euclidian(vector A, vector B){
    // Distance = sqrt(dot(A, A) - 2 * dot(A, B) + dot(B, B))
    return sqrt(arma::sum(arma::pow(A-B, 2)));
}

arma::dvec MatEuclidian(Matrix A, vector B){
    uword col_size = A.n_cols;
    double dist;
    arma::dvec total(col_size, arma::fill::zeros);
    for (uword i = 0; i < col_size; i++) {
        dist = Euclidian(A.col(i), B);
        total(i) = dist;
    }
    return total;
}

arma::dvec MatEuclidian(Matrix A, rvector B){
    uword row_size = A.n_rows;
    double dist;
    arma::dvec total(row_size, arma::fill::zeros);
    for (uword i = 0; i < row_size; i++) {
        dist = Euclidian(A.row(i), B);
        total(i) = dist;
    }
    return total;
}