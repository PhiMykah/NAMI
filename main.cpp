#include <iostream>
#include <string>
#include "Tools/BTS.h"

int main(int argc, char *argv[]) {
    float result;
    float condensed_result;
    int n_atoms = 10;

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

    printf("\nExtended comparison:\n");

    vector ec_results;

    for (int i = 0; i < num_metrics; i++)
    {
        ec_results.push_back(
            ExtendedComparison(matrix, metrics[i], 0, n_atoms)
        );
    }

    for (int i = 0; i < ec_results.size(); i++)
    {
        printf("Metric #%i: %.2f\n", i+1, ec_results[i]);
    }

    printf("\nComplementary Similarity:\n");
    
    Matrix ccs_results;

    for (int i = 0; i < num_metrics; i++){
        ccs_results.push_back(
            CalculateCompSim(matrix, metrics[i], n_atoms)
        );
    }

    for (int i = 0; i < ccs_results.N; i++)
    {
        printf("Metric #%i: ", i+1);
        Matrix(vec2D {ccs_results[i]}).print();
        printf("\n");
    }

    printf("\nCalculate Medoid:\n");

    std::vector<int> cm_results;

    for (int i = 0; i < num_metrics; i++){
        cm_results.push_back(
            CalculateMedoid(matrix, metrics[i], n_atoms)
        );
    }

    for (int i = 0; i < cm_results.size(); i++)
    {
        printf("Metric #%i: ", i+1);
        printf("%i\n", cm_results[i]);
    }

    printf("\nCalculate Outlier:\n");

    std::vector<int> co_results;

    for (int i = 0; i < num_metrics; i++){
        co_results.push_back(
            CalculateOutlier(matrix, metrics[i], n_atoms)
        );
    }

    for (int i = 0; i < co_results.size(); i++)
    {
        printf("Metrix #%i: ", i+1);
        printf("%i\n", co_results[i]);
    }
}