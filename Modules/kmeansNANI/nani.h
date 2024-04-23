#ifndef NANI_H
#define NANI_H 
#include <map>
#include "../../Datatypes/DataContainers.h"

typedef std::map<std::string, float> dict;

struct cluster_data
{   
    cluster_data() {
        labels;
        centers;
        n_iter;
    }
    cluster_data(vector labels_, Matrix centers_, int n_iter_) {
        labels = labels_;
        centers = centers_;
        n_iter = n_iter_;
    }
    vector labels;
    Matrix centers;
    int n_iter;
};

struct scores
{
    scores(){
        ch;
        db;
    }
    scores(float ch_, float db_) {
        ch = ch_;
        db = db_;
    }
    float ch;
    float db;
};

// Initiators for the k-means algorithm
enum class Initiator { COMP_SIM = 0, DIV_SELECT, K_MEANS, RANDOM, VANILLA_KMEANS};

class KmeansNANI
{
public:
    // Constructor
    KmeansNANI(Matrix data, int n_clusters, Metric metric, int n_atoms, Initiator initiator, unsigned short int percentage = 10);

    // Destructor
    ~KmeansNANI();

    void Clear();

    Matrix InitiateKmeans();

    cluster_data KmeansClustering(Matrix initiators);

    cluster_data KmeansClustering(Initiator initiator);

    dict CreateClusterDict(vector labels);

    scores ComputeScores(vector labels);

    void WriteCentroids(Matrix centers, int n_iter);

    cluster_data ExecuteKmeansAll();


private:
    Matrix m_data;
    int n_clusters;
    Metric m_metric;
    int n_atoms;
    Initiator m_initiator;
    int percentage = 10;
    vector m_labels;
    Matrix centers;
    int n_iter;
};

scores ComputeScores(Matrix data, vector labels);
#endif // !NANI_H