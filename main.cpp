#include <iostream>
#include <string>
#include "Tools/BTS.h"

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int n_atoms, float (*test)(Matrix, Metric, int))
{
    printf("%s:\n", title.c_str());
    vector results;

    for (int i = 0; i < metrics.size(); i++) {
        results.push_back(test(matrix, metrics[i], n_atoms));
    }

    for (int i = 0; i < results.size(); i++)
    {
        printf("Metric #%i: %.2f\n", i+1, results[i]);
    }
    printf("\n");
}

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int N, int n_atoms, float (*test)(Matrix, Metric, int, int, float, WFactor))
{
    printf("%s:\n", title.c_str());
    vector results;

    for (int i = 0; i < metrics.size(); i++) {
        results.push_back(test(matrix, metrics[i], N, n_atoms, 0, WFactor::FRACTION));
    }

    for (int i = 0; i < results.size(); i++)
    {
        printf("Metric #%i: %.2f\n", i+1, results[i]);
    }
    printf("\n");
}

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int n_atoms, Matrix (*test)(Matrix, Metric, int))
{
    printf("%s:\n", title.c_str());
    Matrix results;

    for (int i = 0; i < metrics.size(); i++) {
        results.push_back(test(matrix, metrics[i], n_atoms));
    }

    for (int i = 0; i < results.N; i++)
    {
        printf("Metric #%i: ", i+1);
        Matrix(vec2D {results[i]}).print();
        printf("\n");
    }
    printf("\n");
}

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    float percent_trimmed, int n_atoms, Criterion criterion, 
    Matrix (*test)(Matrix, float, Metric, int, Criterion))
{
    printf("%s:\n", title.c_str());
    std::vector<Matrix> results;

    for (int i = 0; i < metrics.size(); i++) {
        results.push_back(test(matrix, percent_trimmed, metrics[i], n_atoms, criterion));
    }

    for (int i = 0; i < results.size(); i++)
    {
        printf("Metric #%i:\n", i+1);
        results[i].print();
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    float result;
    float condensed_result;
    int n_atoms = 10;
    float percent_trimmed = 0.5;
    Criterion criterion = Criterion::COMP_SIM;

    vec2D array {
        { 1.0, 2.0, 3.0, 4.0, 5.0},
        { 6.0, 7.0, 8.0, 9.0,10.0},
        {11.0,12.0,13.0,14.0,15.0},
        {16.0,17.0,18.0,19.0,20.0},
        {21.0,22.0,23.0,24.0,25.0}};

    Matrix matrix(array);

    printf("Matrix:\n");
    matrix.print();
    
    result = MeanSquareDeviation(matrix, n_atoms);

    vector c_sum = matrix.Sum(AXIS::COLUMN);
    vector sq_sum = matrix.pow(2).Sum(AXIS::COLUMN);

    condensed_result = MSDCondensed(c_sum, sq_sum, 5, n_atoms);

    printf("\nMSD Result: %.2f \nMSD Condensed Result: %.2f \n", result, condensed_result);
    
    std::vector<Metric> metrics = {
        Metric::MSD, Metric::BUB, Metric::FAI, 
        Metric::GLE, Metric::JA, Metric::JT, 
        Metric::RT, Metric::RR, Metric::SM, 
        Metric::SS1, Metric::SS2};
    
    int num_metrics = metrics.size();
    
    std::vector<std::string> tests = {
        "Extended Comparison",
        "Complementary Similarity",
        "Calculate Medoid",
        "Calculate Outlier",
        "Trim Outliers"}; 

    // Extended Comparison Results
    OutputResults(tests[0], matrix, metrics, 0, n_atoms, 
    static_cast<float (*)(Matrix, Metric, int, int, float, WFactor)>(&ExtendedComparison));

    // Complementary Similarity Results
    OutputResults(tests[1], matrix, metrics, n_atoms, CalculateCompSim);

    // Calculate Medoid Results
    OutputResults(tests[2], matrix, metrics, n_atoms, 
    [](Matrix mat, Metric met, int atoms)->float{return float(CalculateMedoid(mat, met, atoms));});

    // Calculate Outlier Results
    OutputResults(tests[3], matrix, metrics, n_atoms,
    [](Matrix mat, Metric met, int atoms)->float{return float(CalculateOutlier(mat, met, atoms));});

    // Trim Outliers Results
    OutputResults(tests[4], matrix, metrics, percent_trimmed, n_atoms, criterion,
    static_cast<Matrix (*)(Matrix, float, Metric, int, Criterion)>(&TrimOutliers));
}