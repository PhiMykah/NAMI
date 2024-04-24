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
    char file[] = "./examples/backbone.npy";
    float result;
    float condensed_result;
    int n_atoms = 10;
    float percent_trimmed = 0.5;
    // int percentage = 75;
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
    OutputResults("Extended Comparison", matrix, metrics, 0, n_atoms, 
    static_cast<float (*)(Matrix, Metric, int, int, float, WFactor)>(&ExtendedComparison));

    // Complementary Similarity Results
    OutputResults("Complementary Similarity", matrix, metrics, n_atoms, CalculateCompSim);

    // Calculate Medoid Results
    OutputResults("Calculate Medoid", matrix, metrics, n_atoms, 
    [](Matrix mat, Metric met, int atoms)->float{return float(CalculateMedoid(mat, met, atoms));});

    // Calculate Outlier Results
    OutputResults("Calculate Outlier", matrix, metrics, n_atoms,
    [](Matrix mat, Metric met, int atoms)->float{return float(CalculateOutlier(mat, met, atoms));});

    // Trim Outliers Results
    OutputResults<Criterion>("Trim Outliers", matrix, metrics, percent_trimmed, n_atoms, criterion,
    static_cast<Matrix (*)(Matrix, float, Metric, int, Criterion)>(&TrimOutliers));

    // Trim Outliers Results
    OutputResults<Criterion>("Trim Outliers (MEDOID)", matrix, metrics, percent_trimmed, n_atoms, Criterion::SIM_TO_MEDOID,
    static_cast<Matrix (*)(Matrix, float, Metric, int, Criterion)>(&TrimOutliers));

    // Diversity Selection Medoid Results
    // OutputResults<DiversitySeed>("DiversitySelection + NewIndex (MEDOID)", matrix, metrics,
    // percentage, n_atoms, DiversitySeed::MEDOID, DiversitySelection);

    // Diversity Selection Outlier Results
    // OutputResults<DiversitySeed>("DiversitySelection + NewIndex (OUTLIER)", matrix, metrics,
    // percentage, n_atoms, DiversitySeed::OUTLIER, DiversitySelection);

    // Diversity Selection Random Results
    // OutputResults<DiversitySeed>("DiversitySelection + NewIndex (RANDOM)", matrix, metrics,
    // percentage, n_atoms, DiversitySeed::RANDOM, DiversitySelection);

    // OutputResults<std::vector<int>>("DiversitySelection + NewIndex (LIST)", matrix, metrics,
    // percentage, n_atoms, std::vector<int> {0, 2, 4}, DiversitySelection);
}