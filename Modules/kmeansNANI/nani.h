#ifndef NANI_H
#define NANI_H 
#include "../../Datatypes/DataContainers.h"
#include "../../Tools/BTS/ComplementarySimilarity.h"
#include "../../Tools/BTS/DiversitySelection.h"

typedef arma::field<index_vec> cluster_indices;

struct cluster_data
{   
    cluster_data() {
        labels = vector();
        centers = Matrix();
        n_iter = 0;
    }
    cluster_data(vector labels_, Matrix centers_, uword n_iter_) {
        labels = labels_;
        centers = centers_;
        n_iter = n_iter_;
    }
    vector labels;
    Matrix centers;
    uword n_iter;
};

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

// Initiators for the k-means algorithm
enum class Initiator { COMP_SIM = 0, DIV_SELECT, VANILLA_KMEANS, RANDOM, };

class KmeansNANI
{
public:
    // Constructor
    KmeansNANI(Matrix data, int n_clusters, Metric metric, int n_atoms, Initiator initiator, uword n_iter = 10, unsigned short int percentage = 10);

    // Destructor
    ~KmeansNANI();

    void Clear();

    Matrix InitiateKmeans(Initiator initiator);

    cluster_data KmeansClustering(Matrix initiators);

    cluster_data KmeansClustering();

    cluster_data KmeansClustering(Initiator initiator);

    cluster_indices CreateClusterList(vector labels);

    scores ComputeScores(vector labels);

    void WriteCentroids(Matrix centers);

    // cluster_data ExecuteKmeansAll();

    bool printSteps = false;

private:
    Matrix m_data;
    int n_clusters;
    Metric m_metric;
    int n_atoms;
    Initiator m_initiator;
    int percentage = 10;
    vector m_labels;
    Matrix centers;
    uword n_iter;
};

namespace mlpack {
    class KmeansPlusPlus{
        public:
        //! Empty constructor, required by the InitialPartitionPolicy type definition.
        KmeansPlusPlus();

        /**
         * Initialize the centroids matrix by kmeans++
         * matrix.
         *
         * @param data Dataset.
         * @param clusters Number of clusters.
         * @param centroids Matrix to put initial centroids into.
         */
        template<typename MatType>
        inline static void Cluster(const MatType& data,
                                    const size_t clusters,
                                    arma::mat& centroids)
        {
            centroids.set_size(data.n_rows, clusters);
            // for (size_t i = 0; i < clusters; ++i)
            // {
            //     // Randomly sample a point.
            //     const size_t index = RandInt(0, data.n_cols);
            //     centroids.col(i) = data.col(index);
            // }
        }
    };
}
scores ComputeDataScores(Matrix data, vector labels);
vector GenerateLabels(Matrix data, Matrix centroids);
float CalinskiHarabaszScore(Matrix data, vector labels);
float DaviesBouldinScore(Matrix data, vector labels);
#endif // !NANI_H