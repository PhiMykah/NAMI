#include "../main.h"
#include "../../Modules/kmeansNANI/nani.h"
#include <filesystem>
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    std::string env = getenv("ONE_PIECE");
    std::vector<std::string> confirmations = {"alpha","hairpin","left","pp2"};

    std::vector<std::string> input_files;
    std::string path = env + "/NAMI/examples/ala10";
    for (const auto & entry : fs::directory_iterator(path)) {
        input_files.push_back(entry.path());}

    std::vector<std::string> output_files;
    std::string output_path = env + "/Data/Ala10-tests/NAMI";
    for (const std::string &c : confirmations) {
        output_files.push_back(output_path + "/" + c);}
    

    if (input_files.size() != output_files.size() && input_files.size() != confirmations.size()) {
        std::cout << "Mismatched input and output files!" << std::endl;
        return 1;
    }
        
    Initiator init_types[] = {Initiator::COMP_SIM};
    Metric metric = Metric::MSD;
    int start_n_clusters = 2;
    int end_n_clusters = 30;
    int n_atoms = 109;

    for (size_t i = 0; i != confirmations.size(); i++) 
    {   
        std::cout << "Confirmation: " + confirmations[i] + "\n---------------------" << std::endl;

        // Start of test
        time_t start_time = time(nullptr);

        std::filesystem::path path(output_files[i]);
        std::filesystem::create_directories(path.string());
        std::vector<std::vector<float>> test_results;

        Matrix matrix = loadNPYFile(input_files[i].c_str());

        for (Initiator init_type : init_types)
        {
        
            for (int n_clusters = start_n_clusters; n_clusters < end_n_clusters+1; n_clusters++)
            {
                int n_iter = 20;
                KmeansNANI mod(matrix, n_clusters, metric,
                            n_atoms, init_type, n_iter);

                cluster_data data = mod.KmeansClustering();
                cluster_indices cluster_list = mod.CreateClusterList(data.labels);
                scores results = mod.ComputeScores(data.labels, cluster_list);

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
}