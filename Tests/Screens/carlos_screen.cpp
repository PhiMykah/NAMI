#include "../main.h"
#include "../../Modules/kmeansNANI/nani.h"
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration;

int main(int argc, char *argv[]) {
    char input_file[150];
    char output_dir[150];
    char output[200];
    const char* mat_name = "Data/Carlos/1000_aligned_backbone.npy";
    const char* output_name = "Data/Carlos/nami_outputs";

    sprintf(input_file, "%s/%s", getenv("ONE_PIECE"),mat_name);
    sprintf(output_dir, "%s/%s", getenv("ONE_PIECE"),output_name);
    int n_atoms = 38;

    Initiator init_types[] = {Initiator::COMP_SIM};
    Metric metric = Metric::MSD;
    int start_n_clusters = 2;
    int end_n_clusters = 30;
    int sieve = 1;

    // Start of test
    auto start_time = high_resolution_clock::now();

    std::filesystem::path path(output_dir);
    std::filesystem::create_directories(path.string());
    std::vector<std::vector<float>> test_results;

    Matrix matrix = loadNPYFile(input_file);
    for (Initiator init_type : init_types)
    {
        int n_iter = 20;
        KmeansNANI mod(matrix, start_n_clusters, metric,
                        n_atoms, init_type, n_iter);
        Matrix initial_centroids;// = mod.InitiateKmeans(init_type);
        std::string ic_path = std::string(getenv("ONE_PIECE")) + "/" + std::string("Data/Carlos/initial_centroids.bin");
        initial_centroids.load(ic_path);
        // ASssume n_cols is always greater than or equal to n_clusters
        for (int n_clusters = start_n_clusters; n_clusters < end_n_clusters+1; n_clusters++)
        {
            
            KmeansNANI mod(matrix, n_clusters, metric,
                        n_atoms, init_type, n_iter);

            cluster_data data = mod.KmeansClustering(initial_centroids);
            cluster_indices cluster_list = mod.CreateClusterList(data.labels);
            scores results = mod.ComputeScores(data.centers, cluster_list);

            data.centers.brief_print("Centers:");

            
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

        Matrix output_mat(test_results.at(0).size(), test_results.size());
        for (size_t i = 0; i < test_results.size(); i++)
        {
            output_mat.col(i) = arma::conv_to<vector>::from(test_results.at(i));
        }
        
        output_mat = output_mat.t();

        std::cout << "Init type: " << toStr(init_type)
                  << ", Percentage: " << mod.getPercentage()
                  << ", Metric: " << toStr(metric)
                  << ", Sieve: " << sieve << std::endl;
        
        output_mat.brief_print();

        bool status = true;
        
        sprintf(output, "%s/%d%s_summary.csv", output_dir,mod.getPercentage(),toStr(init_type).c_str());
        arma::field<std::string> header(output_mat.n_cols);
        if (output_mat.n_cols >= 5) {
            header(0) = std::string("Number of clusters") ;
            header(1) = std::string("Number of iterations");
            header(2) = std::string("Calinski-Harabasz score");
            header(3) = std::string("Davies-Bouldin score");
            header(4) = std::string("Average MSD");

            status = output_mat.save(arma::csv_name(output, header));
        } else {
            output_mat.save(arma::csv_name(output, header, arma::csv_opts::no_header));
        }

        if (status == false) {
            std::fprintf(stderr, "Failed to save data to file!\n");
        }

    }

    auto completion_time = high_resolution_clock::now();

    duration<double, std::milli> total_duration = completion_time - start_time;

    std::cout << "Time taken: " << total_duration.count() << std::endl;
}