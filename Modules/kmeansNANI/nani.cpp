#include "nani.h"

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

KmeansNANI::~KmeansNANI()
{
    this->Clear();
}

void KmeansNANI::Clear()
{
    this->m_data = Matrix();
    this->n_clusters = 1;
    this->m_metric = Metric::MSD;
    this->n_atoms = 1;
    this->m_initiator = Initiator::COMP_SIM;
    this->percentage = 100;
}

Matrix KmeansNANI::InitiateKmeans(Initiator initiator)
{
    index_vec initiators_indices;

    if (initiator == Initiator::DIV_SELECT) {
        initiators_indices = DiversitySelection(this->m_data, this->percentage, this->m_metric, DiversitySeed::MEDOID, this->n_atoms);
        return this->m_data.rows(initiators_indices);
    }
    // Kmeans ++ Initialization
    // Comp sim / default
    uword n_total = this->m_data.n_rows;
    uword n_max = (uword)(n_total * this->percentage / 100);
    Matrix comp_sim = CalculateCompSim(this->m_data, this->m_metric, this->n_atoms);
    Matrix sorted_comp_sim = comp_sim;
    if (sorted_comp_sim.n_cols >= 2) {
        // Matrix sorted_comp_sim = sorted(comp_sim, key=lambda item: item[1], reverse=True);
        sortRows(sorted_comp_sim, [](rvector v)->float{return v[1];}, 0, sorted_comp_sim.n_cols-1, true);
    } else {
        fprintf(stderr, "Warning: Comparison Similarity Matrix has size of %llu, \
                        comparison uses first element of vector.", sorted_comp_sim.n_cols);
        sortRows(sorted_comp_sim, [](rvector v)->float{return v[0];}, 0, sorted_comp_sim.n_cols-1, true);
    }

    index_vec total_comp_indices(n_max);


    // vector total_comp_sim_indices = [int(i[0]) for i in sorted_comp_sim][:n_max];
    for (uword i = 0; i < n_max; i++)
    {
        total_comp_indices(i) = (int)sorted_comp_sim(i, 0);
    }

    Matrix top_cc_data = this->m_data.rows(total_comp_indices);
    
    initiators_indices = DiversitySelection(top_cc_data, 100, this->m_metric, DiversitySeed::MEDOID, this->n_atoms);

    return top_cc_data.rows(initiators_indices);
}

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

cluster_data KmeansNANI::KmeansClustering(){
    return this->KmeansClustering(this->m_initiator);
}

cluster_data KmeansNANI::KmeansClustering(Initiator initiator)
{
    Matrix centroids;
    
    switch (initiator)
    {
        case Initiator::DIV_SELECT:
        case Initiator::COMP_SIM:
        case Initiator::VANILLA_KMEANS:
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

cluster_indices KmeansNANI::CreateClusterList(vector labels)
{
    cluster_indices list(this->n_clusters);

    for (int i = 0; i < n_clusters; i++)
    {   
        list(i) = arma::find(labels == i);
    }
    
    return list;
}

scores ComputeDataScores(Matrix data, vector labels)
{
    float ch_score = CalinskiHarabaszScore(data, labels);
    float db_score = DaviesBouldinScore(data, labels);
    return scores(ch_score, db_score);
}

scores KmeansNANI::ComputeScores(vector labels)
{
    return ComputeDataScores(this->m_data, labels);
}

void KmeansNANI::WriteCentroids(Matrix centers)
{
    bool status;
    std::string filename = "centroids.csv";

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

float CalinskiHarabaszScore(Matrix data, vector labels)
{
    return 0.0f;
}

float DaviesBouldinScore(Matrix data, vector labels)
{
    return 0.0f;
}
