#include "main.h"
#include "Modules/kmeansNANI/nani.h"

int main(int argc, char const *argv[])
{
    char file[] = "./examples/backbone.npy";
    int n_clusters = 4;
    int n_atoms = 10;
    uword n_iter = 10;
    ushort percentage = 100;
    Initiator initiator = Initiator::RANDOM;
    Metric metric = Metric::MSD;
    Matrix matrix = loadNPYFile(file);
    KmeansNANI kmn(matrix, n_clusters, metric, 
                   n_atoms, initiator, n_iter, percentage);

    cluster_data data = kmn.KmeansClustering();

    data.labels.t().brief_print("Labels:");
    data.centers.brief_print("Centroids:");
    std::cout << "Max number of iterations: " << data.n_iter << std::endl;

    cluster_indices list = kmn.CreateClusterList(data.labels);

    for (uword i = 0; i < list.n_elem; i++)
    {   
        printf("Cluster #%llu", i);
        list(i).t().brief_print();
    }

    kmn.WriteCentroids(data.centers);
    
    return 0;
}
