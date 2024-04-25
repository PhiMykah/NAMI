#include "nani.h"
#include "../../Tools/BTS/ComplementarySimilarity.h"
#include "../../Tools/BTS/DiversitySelection.h"

KmeansNANI::KmeansNANI(Matrix data, int n_clusters, Metric metric, int n_atoms, Initiator initiator, unsigned short int percentage)
{
    this->m_data = data;
    this->n_clusters = n_clusters;
    this->m_metric = metric;
    this->n_atoms = n_atoms;
    this->m_initiator = initiator;
    this->percentage = percentage;
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

Matrix KmeansNANI::InitiateKmeans()
{
    index_vec initiators_indices;

    switch (this->m_initiator)
    {
    case Initiator::DIV_SELECT:
        initiators_indices = DiversitySelection(this->m_data, this->percentage, this->m_metric, DiversitySeed::MEDOID, this->n_atoms);
        return this->m_data.rows(initiators_indices);
    case Initiator::VANILLA_KMEANS:
        // initiators, indices = kmeans_plusplus(self.data, self.n_clusters, 
        //                                           random_state=None, n_local_trials=1)
        using namespace mlpack;
        arma::mat initiators;
        MatKMeans kmeans;
        kmeans.Cluster(this->m_data, this->n_clusters, initiators);
        
        return arma::conv_to<Matrix>::from(initiators);

    }

    uword n_total = this->m_data.n_rows;
    uword n_max = (uword)(n_total * this->percentage / 100);
    Matrix comp_sim = CalculateCompSim(this->m_data, this->m_metric, this->n_atoms);
    Matrix sorted_comp_sim = comp_sim;
    if (sorted_comp_sim.n_cols >= 2) {
        // Matrix sorted_comp_sim = sorted(comp_sim, key=lambda item: item[1], reverse=True);
        sortRows(sorted_comp_sim, [](rvector v)->float{return v[1];}, 0, sorted_comp_sim.n_cols-1, true);
    } else {
        fprintf(stderr, "Warning: Comparison Similarity Matrix has size of %i, \
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

    return cluster_data();
}

cluster_data KmeansNANI::KmeansClustering(Initiator initiator)
{
    // int n_init = 1;
    // using namespace mlpack;
    // arma::mat initiators;
    // MatKMeans kmeans;
    // kmeans.Cluster(this->m_data, this->n_clusters, initiators);
    return cluster_data();
}

dict KmeansNANI::CreateClusterDict(vector labels)
{
    return dict();
}

scores KmeansNANI::ComputeScores(vector labels)
{
    return scores();
}

void KmeansNANI::WriteCentroids(Matrix centers, int n_iter)
{
}

scores ComputeScores(Matrix data, vector labels)
{
    return scores();
}
