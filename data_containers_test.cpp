#include <iostream>
#include <string>
#include "Datatypes/DataContainers.h"

int main(int argc, char *argv[]) {
    float result;
    float condensed_result;

    vec2D array {
        { 1.0, 2.0, 3.0, 4.0, 5.0},
        { 6.0, 7.0, 8.0, 9.0,10.0},
        {11.0,12.0,13.0,14.0,15.0},
        {16.0,17.0,18.0,19.0,20.0},
        {21.0,22.0,23.0,24.0,25.0}};

    Matrix matrix(array);
    matrix.print();
    std::cout << std::endl;

    Matrix additionA = matrix + 5;
    Matrix additionB = 5 + matrix;
    Matrix subtractionA = matrix - 5;
    Matrix subtractionB = 5 - matrix;
    Matrix mult = matrix * 5;
    Matrix greater = matrix > 12;
    Matrix less = matrix < 12;
    Matrix greater_eq = matrix >= 15.0;
    Matrix less_eq = matrix <= 15.0;

    std::vector<Matrix> matricies {additionA, additionB, subtractionA, subtractionB, mult, greater, less, greater_eq, less_eq};

    for (int i = 0; i < matricies.size(); i++) {
        printf("Matrix #%i\n", i+1);
        matricies[i].print();
        std::cout << std::endl;
    }

    vector sort{20, 3, 2, 13, 4, 5, 9, 19, 10, 16};

    std::vector<int> indicies;
    for (int i = 0; i < sort.size(); i++)
    {
        indicies.push_back(i);
    }

    printf("Unsorted List\n");
    for (int i = 0; i < sort.size(); i++) {
        printf("%.2f ", sort[i]);
    }
    std::cout << std::endl;

    printf("Unsorted indices\n");
    for (int i = 0; i < sort.size(); i++) {
        printf("%i ", indicies[i]);
    }
    std::cout << std::endl;

    quicksort(sort, indicies, 0, sort.size()-1);

    printf("Sorted List\n");
    for (int i = 0; i < sort.size(); i++) {
        printf("%.2f ", sort[i]);
    }
    std::cout << std::endl;

    printf("Sorted indices\n");
    for (int i = 0; i < sort.size(); i++) {
        printf("%i ", indicies[i]);
    }
    std::cout << std::endl;

}