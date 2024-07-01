#ifndef DATA_CONTAINERS_H
#define DATA_CONTAINERS_H
#include <armadillo>
#include <filesystem>
#include <stdexcept>

// Armadillo unsigned integer type
typedef arma::uword uword;

// Armadillo matrix implementation
typedef arma::fmat Matrix;

// Armadillo column vector implementation
typedef arma::fvec vector;

// Armadillo row vector implementation
typedef arma::frowvec rvector;

// Armadillo matrix of indices implementation
typedef arma::uvec index_vec; 

typedef arma::urowvec index_rvec;

// Sum all rows and columns of matrix to singular value
#define MatSum(mat) arma::sum(arma::sum(mat, 0))

// Axis type for 2D matrix
enum class AXIS {COLUMN=0,ROW=1, NONE};

// Metric to use for the extended comparison. Default is 'MSD'.
//     Available metrics:
//     Mean square deviation (MSD), Bhattacharyya's U coefficient (BUB),
//     Faiman's coefficient (Fai), Gleason's coefficient (Gle),
//     Jaccard's coefficient (Ja), Jaccard-Tanimoto coefficient (JT),
//     Rogers-Tanimoto coefficient (RT), Russell-Rao coefficient (RR),
//     Simpson's coefficient (SM), Sokal-Sneath 1 coefficient (SS1),
//     Sokal-Sneath 2 coefficient (SS2).
enum class Metric { MSD=0, BUB, FAI, GLE, JA, JT, RT, RR, SM, SS1, SS2 };

std::string toStr(Metric);

// Type of weight function that will be used. Default is 'FRACTION'.
enum class WFactor { FRACTION=-1, NONE=0};

// Criterion to use for data trimming. Defaults to 'comp_sim'.
// 'comp_sim' removes the most dissimilar objects based on the complement similarity.
// 'sim_to_medoid' removes the most dissimilar objects based on the similarity to the medoid.
enum class Criterion { COMP_SIM=0, SIM_TO_MEDOID };

// Seed of diversity selection. Default is 'MEDOID'.
enum class DiversitySeed { MEDOID=0, OUTLIER, RANDOM, LIST };

// Alignment methods used by functions
enum class AlignMethod { UNI=0, KRON, UNIFORM, KRONECKER };

// ****************************
// Additional Matrix Operations
// ****************************

void sortRows(Matrix mat, float (*key)(rvector v), uword l_index, uword r_index, bool reverse = false); 

// *********************
// Other Data Containers
// *********************

struct scores
{ 
    scores(){
        ch = 0.0;
        db = 0.0;
    }
    scores(float ch_, float db_) {
        ch = ch_;
        db = db_;
    }
    float ch;
    float db;
};

#endif // !DATA_CONTAINERS_H
