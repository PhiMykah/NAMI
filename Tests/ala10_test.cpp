#include "main.h"
#include "../Modules/kmeansNANI/nani.h"

char * outputFile(const char * a, int type) {
    static char output[60];
    switch (type)
    {
    case 1:
        sprintf(output,"output/ala10/NAMI/%s/random_centroids.csv", a);
        break;
    case 2:
        sprintf(output,"output/ala10/NAMI/%s/div_sel_centroids.csv", a);
        break;
    case 0:
    default:
        sprintf(output,"output/ala10/NAMI/%s/comp_sim_centroids.csv", a);
        break;
    }
    
    return output;
}

int main(int argc, char *argv[]) {
    int N_ATOMS = 109;
    int N = 1000;
    int N_CLUSTERS = 4;
    float percent_trimmed = 0.8;
    int percentage = 75;
    int kmn_percentage = 10;
    int n_iter = 33;
    Metric metric = Metric::MSD;

    int n_confirmations = 4;
    const char *data[n_confirmations] = {"./examples/ala10/alpha.npy",
                          "./examples/ala10/hairpin.npy",
                          "./examples/ala10/left.npy",
                          "./examples/ala10/pp2.npy"};
    const char *confirmations[4] = {"alpha", "hairpin", "left", "pp2"};


    int kmeans_test_count = 3;
    std::string tests[kmeans_test_count] = {"******************\nComp Sim Initiator\n******************\n",
                                            "****************\nRandom Initiator\n****************\n",
                                            "********************\nDiv Select Initiator\n********************\n"};

    Initiator initiators[kmeans_test_count] = {Initiator::COMP_SIM, Initiator::RANDOM, Initiator::DIV_SELECT};

    printf("Metric: ");
    std::cout << toStr(metric) << std::endl;
    for (int i = 0; i < n_confirmations; i++)
    {
        Matrix matrix = loadNPYFile(data[i]);
        printf("\nConfirmation: %s\n---------------------\n", confirmations[i]);
        printf("Extended Comparison\n");
        printf("Result: %.4f\n", ExtendedComparison(matrix, metric, N, N_ATOMS, 0, WFactor::FRACTION));

        printf("Complementary Similarity\n");
        CalculateCompSim(matrix, metric, N_ATOMS).t().brief_print("Results:");

        printf("Calculate Medoid\n");
        printf("Result: %d\n", CalculateMedoid(matrix, metric, N_ATOMS));

        printf("Calculate Outlier\n");
        printf("Result: %d\n", CalculateOutlier(matrix, metric, N_ATOMS));
        
        printf("Trim Outliers\n");
        TrimOutliers(matrix, percent_trimmed, metric, N_ATOMS, Criterion::COMP_SIM).brief_print("Results:");
        
        printf("Trim Outliers (MEDOID)\n");
        TrimOutliers(matrix, percent_trimmed, metric, N_ATOMS, Criterion::SIM_TO_MEDOID).brief_print("Results:");
        
        printf("Diversity Selection + NewIndex (MEDOID)\n");
        DiversitySelection(matrix, percentage, metric, DiversitySeed::MEDOID, N_ATOMS).t().brief_print("Results:");

        printf("------\nKMEANS\n------\n");
        for (int j = 0; j < kmeans_test_count; j++)
        {
            std::cout << tests[j];

            KmeansNANI kmn(matrix, N_CLUSTERS, metric, 
                   N_ATOMS, initiators[j], n_iter, kmn_percentage);

            cluster_data data = kmn.KmeansClustering();

            data.labels.t().brief_print("Labels:");

            data.centers.brief_print("Centroids:");
            std::cout << "Max number of iterations: " << data.n_iter << std::endl;

            cluster_indices list = kmn.CreateClusterList(data.labels);

            for (uword p = 0; p < list.n_elem; p++)
            {   
                printf("Cluster #%llu", p);
                list(p).t().brief_print();
            }

            kmn.WriteCentroids(data.centers, outputFile(confirmations[i], j));
        }
        
    }
    
    return 0;
}
