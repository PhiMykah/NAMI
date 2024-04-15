#include "FileIO/ReadNPY.h"
#include "Datatypes/DataContainers.h"

int main(int argc, char const *argv[])
{
    char file[] = "../test.npy";
    Matrix newMat = loadNPYFile(file);
    return 0;
}
