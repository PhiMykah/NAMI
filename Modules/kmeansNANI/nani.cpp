#include "nani.h"
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration;

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

    index_vec total_comp_indices = sorted_comp_sim.subvec(0, n_max-1);

    Matrix top_cc_data = this->m_data.rows(total_comp_indices);
    
    auto t1 = high_resolution_clock::now();
    initiators_indices = DiversitySelection(top_cc_data, 100, this->m_metric, DiversitySeed::MEDOID, this->n_atoms);
    auto t2 = high_resolution_clock::now();
    /* Getting number of milliseconds as a double. */
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "Diversity Selection took " << ms_double.count() << "ms\n";

    Matrix result = top_cc_data.rows(initiators_indices);

    if ((int)result.n_cols < this->n_clusters) {
        throw std::length_error("The number of initiators is less than the number of clusters. Try increasing the percentage.\n");
    }

    return result;
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
cluster_data KmeansNANI::KmeansClustering(Matrix init_centroids)
{
    Matrix centroids = init_centroids.rows(0, this->n_clusters-1).t();
    Matrix data = this->m_data.t();

    arma::kmeans(centroids, data, this->n_clusters, 
            arma::keep_existing, this->n_iter, this->printSteps);

    // Return:
    // - Matrix of this->m_data's shape where every point is the
    //      corresponding group the point is connected to
    // - Centroid matrix
    // - number of max iterations
    vector labels = GenerateLabels(data, centroids);

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
    Matrix data = this->m_data.t();

    switch (initiator)
    {
        case Initiator::DIV_SELECT:
        case Initiator::COMP_SIM:
        case Initiator::VANILLA_KMEANS:
        case Initiator::KMEANS:
            centroids = this->InitiateKmeans(initiator);
            centroids = centroids.rows(0, this->n_clusters-1).t();
            arma::kmeans(centroids, data, this->n_clusters, 
            arma::keep_existing, this->n_iter, this->printSteps);
            break;
        case Initiator::RANDOM:
        default:
            arma::kmeans(centroids, data, this->n_clusters,
            arma::random_subset, this->n_iter, this->printSteps);
            break;
    }

    // Return:
    // - Matrix of this->m_data's shape where every point is the
    //      corresponding group the point is connected to
    // - Centroid matrix
    // - number of max iterations
    vector labels = GenerateLabels(data, centroids);
    
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
        index_vec indicies = arma::find(labels == i);
        list(i) = indicies;
    }
    
    return list;
}

/*
Computes the Davies-Bouldin and Calinski-Harabasz scores.

Parameters
----------
data : Matrix (n_samples, n_features)
    Input dataset.
centers : Matrix
    Matrix of center values (n_features, n_clusters)
clusters : cluster_indices
    Arma field of index vectors corresponding to each cluster
    (n_clusters, variable length)

Returns
-------
scores
    Struct containing the Davies-Bouldin and Calinski-Harabasz scores.
*/
scores ComputeDataScores(Matrix data, Matrix centers, cluster_indices clusters)
{
    float ch_score = CalinskiHarabaszScore(data, centers, clusters);
    float db_score = DaviesBouldinScore(data, centers, clusters);
    return scores(ch_score, db_score);
}

/*
Computes the Davies-Bouldin and Calinski-Harabasz scores using the objects'
internal dataset.

Parameters
----------
centers : Matrix
    Matrix of center values (n_features, n_clusters)
clusters : cluster_indices
    Arma field of index vectors corresponding to each cluster
    (n_clusters, variable length)

Returns
-------
scores
    Struct containing the Davies-Bouldin and Calinski-Harabasz scores.
*/
scores KmeansNANI::ComputeScores(Matrix centers, cluster_indices labels)
{
    return ComputeDataScores(this->m_data, centers, labels);
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

/*
Calculates the Calinski and Harabaz score of the data with given label associations.

Parameters
----------
data : Matrix
    Input dataset.
centers : Matrix
    Matrix of centroid values (n_features, n_clusters)
clusters : cluster_indices
    Arma field of index vectors corresponding to each cluster
    (n_clusters, variable length)

Returns
-------
float
    Calculated Calinski and Harabasz Score.
*/
float CalinskiHarabaszScore(Matrix data, Matrix centers, cluster_indices clusters)
{
    Matrix cluster_k;
    rvector mean = arma::mean(data, COL);
    rvector mean_k;
    float bcss = 0.0f;
    float wcss = 0.0f;

    // Iterate over each cluster and calculate bcss and wcss
    for (uword k = 0; k < centers.n_cols; k++)
    {  
        // Obtain data from clusters
        cluster_k = data.rows(clusters(k));

        // Calculate mean of current data clusters
        mean_k = arma::mean(cluster_k, COL);

        // Calcuate bcss
        bcss += clusters(k).size() * arma::sum(arma::pow((mean_k - mean),2));

        // Calculate wcss
        wcss += arma::accu(arma::pow((cluster_k.each_row() - mean_k),2));
    }
    
    if (wcss == 0) {
        return 1.0;
    }

    return bcss * (data.n_rows - centers.n_cols) / (wcss * (centers.n_cols - 1));
}

/*
Calculates the Davies-Bouldin score of the data with given label associations.

Parameters
----------
data : Matrix
    Input dataset.
labels : vector
    Labels of the k-means algorithm.

Returns
-------
float
    Calculated Davies-Bouldin Score.
*/
float DaviesBouldinScore(Matrix data, Matrix centers, cluster_indices clusters)
{
    /*
    Examples
    --------
    >>> from sklearn.metrics import davies_bouldin_score
    >>> X = [[0, 1], [1, 1], [3, 4]]
    >>> labels = [0, 0, 1]
    >>> davies_bouldin_score(X, labels)
    0.12...
    */

    // Obtain submatrices based on labels
    // X, labels = check_X_y(X, labels)
    // le = LabelEncoder()
    // labels = le.fit_transform(labels)
    
    uword n_features = data.n_cols;

    // n_labels = len(le.classes_)
    uword n_clusters = clusters.n_rows;

    //
    // intra_dists = np.zeros(n_labels)
    std::vector<float>intra_distances(n_clusters);
    
    rvector centroid;
    double average; 
    // centroids = np.zeros((n_labels, len(X[0])), dtype=float)
    Matrix centroids(n_clusters, n_features, arma::fill::zeros);
    // for k in range(n_labels):
    for (uword k = 0; k < n_clusters; k++){
        // cluster_k = _safe_indexing(X, labels == k)
        Matrix data_cluster = data.rows(clusters(k));
        // centroid = cluster_k.mean(axis=0)
        centroid = arma::mean(data_cluster, COL);
        // centroids[k] = centroid
        centroids.row(k) = centroid;

        // intra_dists[k] = np.average(pairwise_distances(cluster_k, [centroid]))
        
        average = arma::mean(MatEuclidian(data_cluster, centroid));
        intra_distances.at(k) = average;
    }

    // centroid_distances = pairwise_distances(centroids)
    // Compare distances between 
    // if np.allclose(intra_dists, 0) or np.allclose(centroid_distances, 0):
        // return 0.0

    // centroid_distances[centroid_distances == 0] = np.inf
    // combined_intra_dists = intra_dists[:, None] + intra_dists
    // scores = np.max(combined_intra_dists / centroid_distances, axis=1)
    // return np.mean(scores)
    return 0.0f;
}

std::string toStr(Initiator init) {
    // enum class Initiator { COMP_SIM = 0, DIV_SELECT, KMEANS, VANILLA_KMEANS, RANDOM, };
    switch (init)
    {
    case Initiator::COMP_SIM:
        return std::string("comp_sim");
    case Initiator::DIV_SELECT:
        return std::string("div_select");
    case Initiator::KMEANS:
        return std::string("kmeans");
    case Initiator::VANILLA_KMEANS:
        return std::string("vanilla_kmeans");
    case Initiator::RANDOM:
        return std::string("random");
    default:
        return std::string("initiator");
    }
}