#include "nani.h"

/*
Default Constructor the KmeansNANI class.

Parameters
----------
data : 2D Matrix (n_samples, n_features)
    Input dataset.
n_clusters : int
    Number of clusters.
metric : Metric enum {'MSD', 'RR', 'JT', etc}
    Metric used for extended comparisons. 
    See `...Datatypes.DataContainers.h` for all available metrics.
initiator : Initiator enum {COMP_SIM, DIV_SELECT, VANILLA_KMEANS, RANDOM}
    'COMP_SIM' selects the inital centers based on the diversity in the densest region of the data.
    'DIV_SELECT' selects the initial centers based on the highest diversity of all data.
    'KMEANS' selects the initial centers based on the greedy k-means++ algorithm.
    'RANDOM' selects the initial centers randomly.
    'VANILLA_KMEANS' selects the initial centers based on the vanilla k-means++ algorithm
n_atoms : int
    Number of atoms. Default is 10.
percentage : int
    Percentage of the dataset to be used for the initial selection of the 
    initial centers. Default is 10.
*/
KmeansNANI::KmeansNANI(Matrix data, int n_clusters, Metric metric, int n_atoms, Initiator initiator, uword n_iter, unsigned short int percentage)
{
    this->m_data = data;
    this->n_clusters = n_clusters;
    this->m_metric = metric;
    this->n_atoms = n_atoms;
    this->m_initiator = initiator;
    this->percentage = percentage;
    this->n_iter = n_iter;
}

/*
KmeansNANI destructor, calls KmeansNANI::Clear()
*/
KmeansNANI::~KmeansNANI()
{
    this->Clear();
}


/*
Reset KmeansNANI object to dummy state
*/
void KmeansNANI::Clear()
{
    this->m_data = Matrix();
    this->n_clusters = 1;
    this->m_metric = Metric::MSD;
    this->n_atoms = 1;
    this->m_initiator = Initiator::COMP_SIM;
    this->percentage = 100;
}


/*
Initializes the k-means algorithm with the selected initiating method
(COMP_SIM, DIV_SELECT, KMEANS, VANILLA_KMEANS).

Defaults to COMP_SIM, RANDOM is handled by the clustering function.

Parameters
----------
initiator : Initiator enum (COMP_SIM, DIV_SELECT, KMEANS, VANILLA_KMEANS)


Returns
-------
Matrix
    The initial centers for k-means of shape (n_features, n_clusters).
*/
Matrix KmeansNANI::InitiateKmeans(Initiator initiator)
{
    index_vec initiators_indices;

    if (initiator == Initiator::DIV_SELECT) {
        initiators_indices = DiversitySelection(this->m_data, this->percentage, this->m_metric, DiversitySeed::MEDOID, this->n_atoms);
        return this->m_data.rows(initiators_indices);
    }
    // Kmeans++ Initialization
    // Vanilla Kmeans++ Initialization
    // Comp sim / default
    uword n_total = this->m_data.n_rows;
    uword n_max = (uword)(n_total * this->percentage / 100);
    vector comp_sim = CalculateCompSim(this->m_data, this->m_metric, this->n_atoms);
    index_vec sorted_comp_sim = arma::sort_index(comp_sim, "descend");

    index_vec total_comp_indices(n_max);

    // vector total_comp_sim_indices = [int(i[0]) for i in sorted_comp_sim][:n_max];
    for (uword i = 0; i < n_max; i++)
    {
        total_comp_indices(i) = sorted_comp_sim(i);
    }

    Matrix top_cc_data = this->m_data.rows(total_comp_indices);
    
    initiators_indices = DiversitySelection(top_cc_data, 100, this->m_metric, DiversitySeed::MEDOID, this->n_atoms);

    return top_cc_data.rows(initiators_indices);
}


/*
Executes the k-means algorithm with given initial centroid matrix.

Parameters
----------
initiators : 2D Matrix of (n_features, n_clusters)    

Returns
-------
cluster_data
    Struct containing:
        - The labels of each point to the closest centroid
        - Matrix of centroids
        - Maximum allowed iterations.
*/
cluster_data KmeansNANI::KmeansClustering(Matrix initiators)
{
    
    Matrix centroids(initiators);

    arma::kmeans(centroids, this->m_data, this->n_clusters,
                 arma::keep_existing, this->n_iter, this->printSteps);


    // Return:
    // - Matrix of this->m_data's shape where every point is the
    //      corresponding group the point is connected to
    // - Centroid matrix
    // - number of max iterations
    vector labels = GenerateLabels(this->m_data, centroids);

    return cluster_data(labels, centroids, this->n_iter);
}


/*
Executes KmeansClustering(Initiator initiator) with initiator stored by class

Returns
-------
cluster_data
    Struct containing:
        - The labels of each point to the closest centroid
        - Matrix of centroids
        - Maximum allowed iterations.
*/
cluster_data KmeansNANI::KmeansClustering(){
    return this->KmeansClustering(this->m_initiator);
}


/*
Executes the k-means algorithm with given initiator to generate
starting centroid

Parameters
----------
initiators : Initiator enum (COMP_SIM, DIV_SELECT, KMEANS, VANILLA_KMEANS, RANDOM)

Returns
-------
cluster_data
    Struct containing:
        - The labels of each point to the closest centroid
        - Matrix of centroids
        - Maximum allowed iterations.
*/
cluster_data KmeansNANI::KmeansClustering(Initiator initiator)
{
    Matrix centroids;
    
    switch (initiator)
    {
        case Initiator::DIV_SELECT:
        case Initiator::COMP_SIM:
        case Initiator::VANILLA_KMEANS:
        case Initiator::KMEANS:
            centroids = this->InitiateKmeans(initiator);
            arma::kmeans(centroids, this->m_data, this->n_clusters, 
            arma::keep_existing, this->n_iter, this->printSteps);
            break;
        case Initiator::RANDOM:
        default:
            arma::kmeans(centroids, this->m_data, this->n_clusters,
            arma::random_subset, this->n_iter, this->printSteps);
            break;
    }

    // Return:
    // - Matrix of this->m_data's shape where every point is the
    //      corresponding group the point is connected to
    // - Centroid matrix
    // - number of max iterations
    vector labels = GenerateLabels(this->m_data, centroids);
    
    return cluster_data(labels, centroids, this->n_iter);
}

/*
Creates a mapping between the cluster labels and
the vector of indices that correspond to the cluster.

Parameters
----------
labels : vector
    Labels of the k-means algorithm.

Returns
-------
cluster_indices
    arma::field map with labels as keys and the indices of the data as values.
*/
cluster_indices KmeansNANI::CreateClusterList(vector labels)
{
    cluster_indices list(this->n_clusters);

    for (int i = 0; i < n_clusters; i++)
    {   
        list(i) = arma::find(labels == i);
    }
    
    return list;
}

/*
Computes the Davies-Bouldin and Calinski-Harabasz scores.

Parameters
----------
data : Matrix (n_samples, n_features)
    Input dataset.
labels : vector
    Labels of the k-means algorithm.

Returns
-------
scores
    Struct containing the Davies-Bouldin and Calinski-Harabasz scores.
*/
scores ComputeDataScores(Matrix data, vector labels)
{
    float ch_score = CalinskiHarabaszScore(data, labels);
    float db_score = DaviesBouldinScore(data, labels);
    return scores(ch_score, db_score);
}

/*
Computes the Davies-Bouldin and Calinski-Harabasz scores using the objects'
internal dataset.

Parameters
----------
labels : vector
    Labels of the k-means algorithm.

Returns
-------
scores
    Struct containing the Davies-Bouldin and Calinski-Harabasz scores.
*/
scores KmeansNANI::ComputeScores(vector labels)
{
    return ComputeDataScores(this->m_data, labels);
}

/*
Writes the centroids of the k-means algorithm to a file.

Parameters
----------
centers : 2D Matrix (n_features, n_clusters)
    Centroids of the k-means algorithm.
filename : std::string 
    String representation of output file path (including extension).
*/
void KmeansNANI::WriteCentroids(Matrix centers, std::string filename)
{
    if (centers.is_empty()) {return;}
    
    bool status = true;

    std::filesystem::path path(filename);

    if (path.has_parent_path()) {
        status = std::filesystem::create_directories(path.parent_path().string());
    }

    arma::field<std::string> header(centers.n_cols);
    if (centers.n_cols >= 2) {
    header(0) = std::string("Number of clusters: ") + std::to_string(this->n_clusters);
    header(1) = std::string("Number of iterations: ") + std::to_string(this->n_iter);
    status = centers.save(arma::csv_name(filename, header));
    } else {
        centers.save(arma::csv_name(filename, header, arma::csv_opts::no_header));
    }

    if (status == false) {
        std::fprintf(stderr, "Failed to save centroids to file!\n");
    }
}

/*
Generate vector of centroid labels based on the closest centroid to each data point

Parameters
----------
data : 2D Matrix (n_samples, n_features)
    Input dataset
centroids : 2D Matrix (n_features, n_clusters)
    center values 
Returns
-------
vector
    vector of center labels corresponding to each sample
*/
vector GenerateLabels(Matrix data, Matrix centroids)
{
    // Create a label list with number of features
    vector labelVector(data.n_cols);

    // Initialize variables 
    float minimum_distance; // Minimal distance for label calculation
    vector square_diff; // Squared difference between all elements 
    float distance; // Distance calculated for the given center
    int label; // Corresponding label index for element

    // Iterate over all features and their corresponding labelVector index
    for (uword i = 0; i < data.n_cols; i++)
    {
        // Initialize minimum distance 
        minimum_distance = MAXFLOAT;
        for (uword center = 0; center < centroids.n_cols; center++)
        {
            // Calculate euclidian distance between center and feature
            square_diff = arma::pow((centroids.col(center) - data.col(i)), 2);
            distance = std::sqrt(arma::sum(square_diff));

            // Set label to current center index if the distance is less than minimum
            if (distance < minimum_distance) {
                minimum_distance = distance;
                label = center;
            }

        }
        labelVector(i) = label;
    }

    return labelVector;
}