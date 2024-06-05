#include "main.h"
#include "../Modules/kmeansNANI/nani.h"

int main(int argc, char const *argv[])
{
    char file[] = "../examples/backbone.npy";
    int n_clusters = 4;
    int n_atoms = 10;
    uword n_iter = 10;
    ushort percentage = 10;
    int test_count = 3;
    std::string tests[test_count] = {"****************\nRandom Initiator\n****************\n",
                                    "******************\nComp Sim Initiator\n******************\n",
                                   "********************\nDiv Select Initiator\n********************\n"};
    std::string rc = "output/random_centroids.csv";
    std::string csc = "output/comp_sim_centroids.csv";
    std::string dsc = "output/div_sel_centroids.csv";
    std::string filenames[test_count] = {rc, csc, dsc};

    Initiator initiators[test_count] = {Initiator::RANDOM, Initiator::COMP_SIM, Initiator::DIV_SELECT};
    Metric metric = Metric::MSD;
    Matrix matrix = loadNPYFile(file);

    for (int i = 0; i < test_count; i++)
    {   
        std::cout << tests[i];

        KmeansNANI kmn(matrix, n_clusters, metric, 
                   n_atoms, initiators[i], n_iter, percentage);

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

        kmn.WriteCentroids(data.centers, filenames[i]);
    }

    return 0;
}
