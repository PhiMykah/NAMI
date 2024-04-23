#include "main.h"

int main(int argc, char const *argv[])
{
    char file[] = "./examples/backbone.npy";
    int n_atoms = 10;

    Matrix matrix = loadNPYFile(file);
    float result = ExtendedComparison(matrix, Metric::MSD, 0, n_atoms, 0, WFactor::FRACTION);

    std::cout << result << std::endl;
    return 0;
}
