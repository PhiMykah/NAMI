#include "main.h"

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int n_atoms, float (*test)(Matrix, Metric, int))
{
    printf("%s:\n", title.c_str());
    rvector results(metrics.size());

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        results(i) = test(matrix, metrics[i], n_atoms);
    }
    
    for (long unsigned int i = 0; i < results.size(); i++)
    {
        printf("Metric \'%s\': %.2f\n", toStr(metrics[i]).c_str(), results[i]);
    }
    printf("\n");
}

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int N, int n_atoms, float (*test)(Matrix, Metric, int, int, float, WFactor))
{
    printf("%s:\n", title.c_str());
    rvector results(metrics.size());

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        results(i) = test(matrix, metrics[i], N, n_atoms, 0, WFactor::FRACTION);
    }

    for (long unsigned int i = 0; i < results.size(); i++)
    {
        printf("Metric \'%s\': %.2f\n", toStr(metrics[i]).c_str(), results[i]);
    }
    printf("\n");
}

void OutputResults( 
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int n_atoms, Matrix (*test)(Matrix, Metric, int))
{
    printf("%s:\n", title.c_str());
    std::vector<Matrix> results;

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        results.push_back(test(matrix, metrics[i], n_atoms).as_row());
    }

    for (long unsigned int i = 0; i < metrics.size(); i++)
    {
        printf("Metric \'%s\':\n", toStr(metrics[i]).c_str());
        results[i].brief_print();
        printf("\n");
    }
    printf("\n");
}

template <typename Enum> void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    float percent_trimmed, int n_atoms, Enum criterion, 
    Matrix (*test)(Matrix, float, Metric, int, Enum))
{
    printf("%s:\n", title.c_str());
    std::vector<Matrix> results;

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        results.push_back(test(matrix, percent_trimmed, metrics[i], n_atoms, criterion));
    }

    for (long unsigned int i = 0; i < results.size(); i++)
    {
        printf("Metric \'%s\':\n", toStr(metrics[i]).c_str());
        results[i].brief_print();
        printf("\n");
    }
    printf("\n");
}

template <typename T> void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int percentage, int n_atoms, T start,
    std::vector<int> (*test)(Matrix, int, Metric, T, int))
{
    printf("%s:\n", title.c_str());
    std::vector<std::vector<int>> results;
    
    for (long unsigned int i = 0; i < metrics.size(); i++) {
        results.push_back(test(matrix, percentage, metrics[i], start, n_atoms));
    }
        for (long unsigned int i = 0; i < results.size(); i++)
    {
        printf("Metric \'%s\':\n", toStr(metrics[i]).c_str());
        
        for (long unsigned int j = 0; j < results[i].size(); j++)
        {
            printf("%i ", results[i][j]);
        }
        printf("\n\n");
    }
}

int main(int argc, char *argv[]) {
    char file[] = "../examples/backbone.npy";
    float result;
    float condensed_result;
    int n_atoms = 10;
    float percent_trimmed = 0.5;
    int percentage = 75;
    Criterion criterion = Criterion::COMP_SIM;

    Matrix matrix = loadNPYFile(file);

    // printf("Matrix:\n");
    // matrix.print();
    
    result = MeanSquareDeviation(matrix, n_atoms);

    rvector c_sum = arma::sum(matrix,(int)AXIS::COLUMN);
    rvector sq_sum = arma::sum(arma::pow(matrix,2),(int)AXIS::COLUMN);

    condensed_result = MSDCondensed(c_sum, sq_sum, matrix.n_rows, n_atoms);

    printf("\nMSD Result: %.2f \nMSD Condensed Result: %.2f \n", result, condensed_result);
    
    std::vector<Metric> metrics = {
        Metric::MSD, Metric::BUB, Metric::FAI, 
        Metric::GLE, Metric::JA, Metric::JT, 
        Metric::RT, Metric::RR, Metric::SM, 
        Metric::SS1, Metric::SS2};
    
    // Extended Comparison Results
    printf("Extended Comparison\n");
    rvector ec_results(metrics.size());

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        ec_results(i) = ExtendedComparison(matrix, metrics[i], 0, n_atoms, 0, WFactor::FRACTION);
    }

    for (long unsigned int i = 0; i < ec_results.size(); i++)
    {
        printf("Metric \'%s\': %.2f\n", toStr(metrics[i]).c_str(), ec_results[i]);
    }
    printf("\n");
    
    
    // Complementary Similarity Results
    printf("Complementary Similarity\n");
    std::vector<Matrix> cs_results;

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        cs_results.push_back(CalculateCompSim(matrix, metrics[i], n_atoms).as_row());
    }

    for (long unsigned int i = 0; i < metrics.size(); i++)
    {
        printf("Metric \'%s\':\n", toStr(metrics[i]).c_str());
        cs_results[i].brief_print();
        printf("\n");
    }
    printf("\n");

    // Calculate Medoid Results
    printf("Calculate Medoid\n");
    rvector cm_results(metrics.size());

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        cm_results(i) = CalculateMedoid(matrix, metrics[i], n_atoms);
    }
    
    for (long unsigned int i = 0; i < cm_results.size(); i++)
    {
        printf("Metric \'%s\': %.2f\n", toStr(metrics[i]).c_str(), cm_results[i]);
    }
    printf("\n");

    // Calculate Outlier Results
    printf("Calculate Outlier\n");
    rvector co_results(metrics.size());

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        co_results(i) = CalculateOutlier(matrix, metrics[i], n_atoms);
    }
    
    for (long unsigned int i = 0; i < co_results.size(); i++)
    {
        printf("Metric \'%s\': %.2f\n", toStr(metrics[i]).c_str(), co_results[i]);
    }
    printf("\n");

    // Trim Outliers Results
    printf("Trim Outliers\n");
    std::vector<Matrix> to_results;

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        to_results.push_back(TrimOutliers(matrix, percent_trimmed, metrics[i], n_atoms, criterion));
    }

    for (long unsigned int i = 0; i < to_results.size(); i++)
    {
        printf("Metric \'%s\':\n", toStr(metrics[i]).c_str());
        to_results[i].brief_print();
        printf("\n");
    }
    printf("\n");

    // Trim Outliers Results
    printf("Trim Outliers (MEDOID)\n");
    std::vector<Matrix> tom_results;

    for (long unsigned int i = 0; i < metrics.size(); i++) {
        tom_results.push_back(TrimOutliers(matrix, percent_trimmed, metrics[i], n_atoms, Criterion::SIM_TO_MEDOID));
    }

    for (long unsigned int i = 0; i < tom_results.size(); i++)
    {
        printf("Metric \'%s\':\n", toStr(metrics[i]).c_str());
        tom_results[i].brief_print();
        printf("\n");
    }
    printf("\n");


    // Diversity Selection Medoid Results
    printf("Diversity Selection + NewIndex (MEDOID)\n");
    std::vector<index_vec> ds_results;
    
    percentage = 10;
    for (long unsigned int i = 0; i < metrics.size(); i++) {
        ds_results.push_back(DiversitySelection(matrix, percentage, metrics[i], DiversitySeed::MEDOID, n_atoms));
    }

    for (long unsigned int i = 0; i < ds_results.size(); i++) {
        printf("Metric \'%s\':\n", toStr(metrics[i]).c_str());
        ds_results[i].t().brief_print();
        printf("\n");
    }
    printf("\n");


    // Diversity Selection Outlier Results
    // OutputResults<DiversitySeed>("DiversitySelection + NewIndex (OUTLIER)", matrix, metrics,
    // percentage, n_atoms, DiversitySeed::OUTLIER, DiversitySelection);

    // Diversity Selection Random Results
    // OutputResults<DiversitySeed>("DiversitySelection + NewIndex (RANDOM)", matrix, metrics,
    // percentage, n_atoms, DiversitySeed::RANDOM, DiversitySelection);

    // OutputResults<std::vector<int>>("DiversitySelection + NewIndex (LIST)", matrix, metrics,
    // percentage, n_atoms, std::vector<int> {0, 2, 4}, DiversitySelection);
}
