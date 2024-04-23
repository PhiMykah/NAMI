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
}

Matrix KmeansNANI::InitiateKmeans()
{
    std::vector<int> initiators_indices;

    switch (this->m_initiator)
    {
    case Initiator::DIV_SELECT:
        initiators_indices = DiversitySelection(this->m_data, this->percentage, this->m_metric, DiversitySeed::MEDOID, this->n_atoms);
        return this->m_data[initiators_indices];
    case Initiator::VANILLA_KMEANS:
        // initiators, indices = kmeans_plusplus(self.data, self.n_clusters, 
        //                                           random_state=None, n_local_trials=1)
        break;
    }

    int n_total = this->m_data.M;
    int n_max = (int)(n_total * this->percentage / 100);
    Matrix comp_sim = CalculateCompSim(this->m_data, this->m_metric, this->n_atoms);
    Matrix sorted_comp_sim = comp_sim;
    if (sorted_comp_sim.M >= 2) {
        // Matrix sorted_comp_sim = sorted(comp_sim, key=lambda item: item[1], reverse=True);
        sorted(sorted_comp_sim, [](vector v)->float{return v[1];}, 0, sorted_comp_sim.M-1, true);
    } else {
        fprintf(stderr, "Warning: Comparison Similarity Matrix has size of %i, \
                        comparison uses first element of vector.", sorted_comp_sim.M);
        sorted(sorted_comp_sim, [](vector v)->float{return v[0];}, 0, sorted_comp_sim.M-1, true);
    }

    std::vector<int> total_comp_indices;
    int i = 0;

    // vector total_comp_sim_indices = [int(i[0]) for i in sorted_comp_sim][:n_max];
    while (i < n_max)
    {
        total_comp_indices.push_back((int) sorted_comp_sim[i][0]);
    }

    Matrix top_cc_data = this->m_data[total_comp_indices];
    
    initiators_indices = DiversitySelection(top_cc_data, 100, this->m_metric, DiversitySeed::MEDOID, this->n_atoms);

    return top_cc_data[initiators_indices];
}

cluster_data KmeansNANI::KmeansClustering(Matrix initiators)
{
    return cluster_data();
}

cluster_data KmeansNANI::KmeansClustering(Initiator initiator)
{
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
