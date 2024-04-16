#include "FileIO/ReadNPY.h"
#include "Datatypes/DataContainers.h"

int main(int argc, char const *argv[])
{   
    char file[] = "../backbone.npy";
    Matrix newMat = loadNPYFile(file);
    newMat.print();
    return 0;
}
