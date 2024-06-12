#include "main.h"
#include "../Modules/kmeansNANI/nani.h"

int main(int argc, char *argv[]) {
    char input_file[150];
    char output_dir[150];
    const char* mat_name = "Data/Carlos/1000_aligned_backbone.npy";
    const char* output_name = "Data/Carlos/nami_outputs";

    sprintf(input_file, "%s/%s", getenv("ONE_PIECE"),mat_name);
    sprintf(output_dir, "%s/%s", getenv("ONE_PIECE"),output_name);
    int n_atoms = 38;

    Initiator init_types[] = {Initiator::COMP_SIM};
    Metric metric = Metric::MSD;
    int start_n_clusters = 2;
    int end_n_clusters = 30;

    // Start of test
    time_t start_time = time(nullptr);

    std::filesystem::path path(output_dir);
    std::filesystem::create_directories(path.string());
    std::vector<std::vector<float>> test_results;

    Matrix matrix = loadNPYFile(input_file);
    for (Initiator init_type : init_types)
    {
        
        for (int n_clusters = start_n_clusters; n_clusters < end_n_clusters; n_clusters++)
        {
            int n_iter = 20;
            KmeansNANI mod(matrix, n_clusters, metric,
                        n_atoms, init_type, n_iter);

            cluster_data data = mod.KmeansClustering();
            scores results = mod.ComputeScores(data.labels);

            cluster_indices cluster_list = mod.CreateClusterList(data.labels);
            
            float msd_total = 0;
            float msd;
            Matrix sub_mat;
            for (uword i = 0; i < cluster_list.size(); i++)
            {
                sub_mat = matrix.rows(cluster_list[i]);
                msd = ExtendedComparison(sub_mat, metric, 0, n_atoms);

                msd_total += msd;
            }
            std::vector<float> all_scores{float(cluster_list.size()), float(n_iter), 
                                  results.ch, results.db, msd_total/cluster_list.size()};
            test_results.push_back(all_scores);
        }

    }

    // Matrix output_mat = arma::conv_to< Matrix >::from(test_results);
    // output_mat.brief_print();

    time_t completion_time = time(nullptr);

    double duration = difftime(completion_time, start_time);

    std::cout << "Time taken: " << duration << std::endl;
}