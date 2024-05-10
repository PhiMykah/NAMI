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

// Initiators for the k-means algorithm
enum class Initiator { COMP_SIM = 0, DIV_SELECT, KMEANS, VANILLA_KMEANS, RANDOM, };

/*
K-means algorithm with the N-Ary Natural Initialization (NANI).
    
Attributes
----------
m_data : 2D matrix (n_samples, n_features)
    Input dataset.
n_clusters : int
    Number of clusters.
m_metric : Metric enum {'MSD', 'RR', 'JT', etc}
    Metric used for extended comparisons. 
    See `...Datatypes.DataContainers.h` for all available metrics.
n_atoms : int
    Number of atoms.
m_initiator : Initiator enum {COMP_SIM, DIV_SELECT, VANILLA_KMEANS, RANDOM}
    Type of initiator selection. 
n_iter : int
    Max number of iterations run.
percentage : int
    Percentage of the dataset to be used for the initial selection of the 
    initial centers. Default is 10.
m_labels : vector of length n_samples
    Labels of each point.
centers : 2D matrix (n_clusters, n_features)
    Cluster centers.
cluster_dict : dict
    Dictionary of the clusters and their corresponding indices.
*/
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

    void WriteCentroids(Matrix centers, std::string filename = "centroids.csv");

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
#endif // !NANI_H