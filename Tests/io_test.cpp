#include "../FileIO/ReadNPY.h"
#include "../Datatypes/DataContainers.h"

int main(int argc, char const *argv[])
{   
    char file[] = "../examples/backbone.npy";
    Matrix newMat = loadNPYFile(file);
    // std::cout << newMat.n_rows << " " << newMat.n_cols << std::endl;
    newMat.print();
    return 0;
}
