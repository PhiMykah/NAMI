#ifndef READ_NPY_H
#define READ_NPY_H
#include <fstream>
#include <iostream>
#include <string>
#include "../Datatypes/DataContainers.h"

enum class DataType {f8=0, f4=1, i8=2, i4=3};

// Double size
#define D_SIZE sizeof(double)

// Float size
#define F_SIZE sizeof(float)

// int32 size
#define I4_SIZE sizeof(int32_t)

// int64 size
#define I8_SIZE sizeof(int64_t)

std::string toString(char* a, int size);

struct HeaderNPY
{
    HeaderNPY(
        DataType dtype_,
        bool fortran_,
        int M_size_, int N_size_)
    {
        this->dtype = dtype_;
        this->fortran_order = fortran_;
        M_size = M_size_;
        N_size = N_size_;
    }

    DataType dtype;
    bool fortran_order;
    int M_size;
    int N_size;
};

Matrix loadNPYFile(char file_path[]);

HeaderNPY parseHeader(std::string header, int header_size);

Matrix PopulateMatrix(std::ifstream &array_file, HeaderNPY header, std::streampos data_start);

template <typename DType> DType GetDTypeFromBytes(char value[sizeof(DType)]);

#endif // !READ_NPY_H