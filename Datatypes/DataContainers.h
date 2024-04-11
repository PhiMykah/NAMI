#ifndef DATA_CONTAINERS_H
#define DATA_CONTAINERS_H
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <limits.h>
#include <cstdlib>

// float vector definition
typedef std::vector<float> vector;

// Vector of float vectors (2D vector)
typedef std::vector<vector> vec2D;

// Empty vector object
static const vector EMPTY_VECTOR;

// Axis type for 2D matrix
enum class AXIS {COLUMN=0,ROW=1, NONE};

// Metric to use for the extended comparison. Default is 'MSD'.
//     Available metrics:
//     Mean square deviation (MSD), Bhattacharyya's U coefficient (BUB),
//     Faiman's coefficient (Fai), Gleason's coefficient (Gle),
//     Jaccard's coefficient (Ja), Jaccard-Tanimoto coefficient (JT),
//     Rogers-Tanimoto coefficient (RT), Russell-Rao coefficient (RR),
//     Simpson's coefficient (SM), Sokal-Sneath 1 coefficient (SS1),
//     Sokal-Sneath 2 coefficient (SS2).
enum class Metric { MSD=0, BUB, FAI, GLE, JA, JT, RT, RR, SM, SS1, SS2 };

// Type of weight function that will be used. Default is 'FRACTION'.
enum class WFactor { FRACTION=-1, NONE=0};

// Criterion to use for data trimming. Defaults to 'comp_sim'.
// 'comp_sim' removes the most dissimilar objects based on the complement similarity.
// 'sim_to_medoid' removes the most dissimilar objects based on the similarity to the medoid.
enum class Criterion { COMP_SIM=0, SIM_TO_MEDOID };

// Seed of diversity selection. Default is 'MEDOID'.
enum class DiversitySeed { MEDOID=0, OUTLIER, RANDOM, LIST };

// Alignment methods used by functions
enum class AlignMethod { UNI=0, KRON, UNIFORM, KRONECKER };

// Quicksort helper class
void quicksort(vector &vec, int l_index, int r_index);

// Quicksort that also sorts an index list
void quicksort(vector &vec, std::vector<int> &indx_vector, int l_index, int r_index);

// 2D-Matrix object
class Matrix
{
public:
    // Default Constructor
    Matrix();

    // Overloaded Constructor
    Matrix(vec2D array);

    // Overloaded Constructor with existent flat
    Matrix(vec2D array, vector flat);

    // Array member setter
    void SetArray(vec2D m_array);

    // Array member getter
    vec2D GetArray();
    
    // Remove vector at index
    void erase(int);

    // Remove vectors at indicies
    void erase(std::vector<int>);

    void vSwap(int l_index, int r_index);
    // ********************
    // Overloaded Operators
    // ********************

    // Conversion operators
    // ====================

    operator vector() const;

    // Index matrix rows
    vector const& operator[](int) const;
    
    // Return all matching values from indicies matrix
    Matrix operator[](Matrix) const;

    // Arithmetic 
    // ==========

    // Addition
    Matrix operator+(Matrix x);

    // Subtraction
    Matrix operator-(Matrix x);

    // Multiplication
    Matrix operator*(Matrix x);

    // Division
    Matrix operator/(Matrix x);

    // Modulo operation
    Matrix operator%(Matrix x);


    // Elementwise Arithmetic 
    // ======================

    // Addition elementwise
    Matrix operator+(float x);

    // Subtraction elementwise
    Matrix operator-(float x);

    // Multiplication elementwise 
    Matrix operator* (float x);

    // Division elementwise
    Matrix operator/ (float x);

    // Modulo operation elementwise
    Matrix operator% (float x);

    // Vectorwise Arithmetic 
    // =====================

    // Addition vectorwise
    Matrix operator+(vector x);

    // Subtraction vectorwise
    Matrix operator-(vector x);

    // Multiplication vectorwise
    Matrix operator*(vector x);

    // Division vectorwise
    Matrix operator/(vector x);

    // Modulo operation elementwise
    Matrix operator%(vector x);

    // **************** 
    // vec2D Operations 
    // ****************

    // Return the absolute value matrix
    Matrix Absolute();

    // Flatten the 2D array into a 1D array and add to flatten member
    void Flatten();

    // Summate all of the values overall
    float Sum();

    // Summate all of the values along a particular axis
    vector Sum(AXIS axis);

    // Perform an elementwise operation on the data
    Matrix Elementwise(float (*func)(float)) const;

    // Perform an elementwise operation on the data with float input
    Matrix Elementwise(Matrix x, float (*func)(float, float)) const;

    // Add another vector to the matrix
    void push_back(vector x);

    // Raise each element of the matrix to the input power
    Matrix pow(int power);

    // Temporary print function until implemented through outstream
    void print();

    // Called by destructor to clear object
    void clear();

    // Destructor
    ~Matrix();

    // Number of rows
    int N;  
    // Number of columns
    int M;

private:
    // Container where the data is stored
    vec2D m_array; 
    vector m_flat;
};

void quicksortMatrix(Matrix &mat, std::vector<int> &scores, int l_index, int r_index);

// *************************
// * Reshaping and Casting *
// *************************

// Cast Constant value to Matrix
Matrix cast(float x, const Matrix& source_matrix); 

// Cast vector of floats to a 2D vector
Matrix cast(vector vec, const Matrix& source_matrix);

// Reshape Vector to Matrix
Matrix vectorReshape(vector vec, const Matrix& source_matrix);

// Function Lambdas for elementwise operations
#define ADD [](float a, float b)->float{return a + b;}
#define SUB [](float a, float b)->float{return a - b;}
#define MULT [](float a, float b)->float{return a * b;}
#define DIV [](float a, float b)->float{return a / b;}
#define LESS_THAN [](float a, float b)->float{if (a < b) {return 1.0;} return 0.0;}
#define GREATER_THAN [](float a, float b)->float{if (a > b) {return 1.0;} return 0.0;}
#define LEQ [](float a, float b)->float{if (a <= b) {return 1.0;} return 0.0;}
#define GEQ [](float a, float b)->float{if (a >= b) {return 1.0;} return 0.0;}
#define POW [](float a, float b)->float{return std::pow(a,b);}
#define MOD [](float a, float b)->float{return (int) a % (int) b;}

// *************************
// More Overloaded Operators
// *************************

// Matrix Arithmetic
Matrix pow(const Matrix& x, const Matrix& y);
Matrix pow(const Matrix& x, float y);
Matrix pow(const Matrix& x, vector y);

// Vector-Vector arithmetic
vector operator+(vector x, vector y);
vector operator-(vector x, vector y);
vector pow(vector x, float y);

// Matrix-Constant arithmetic
Matrix operator+(float x, const Matrix& y);
Matrix operator-(float x, const Matrix& y);
Matrix operator* (float x, const Matrix& y);
Matrix operator/ (float x, const Matrix& y);
Matrix pow(float x, const Matrix& y);

// Matrix-Constant inequalities
Matrix operator <(const Matrix& lhs, const float rhs);
Matrix operator >(const Matrix& lhs, const float rhs);
Matrix operator <=(const Matrix& lhs, const float rhs);
Matrix operator >=(const Matrix& lhs, const float rhs);

// Matrix-Vector arithmetic
Matrix operator+(vector x, const Matrix& y);
Matrix operator-(vector x, const Matrix& y);
Matrix operator* (vector x, const Matrix& y);
Matrix operator/ (vector x, const Matrix& y);
Matrix pow(vector x, const Matrix& y);

// Matrix-Vector inequalities
Matrix operator <(const Matrix& lhs, const vector rhs);
Matrix operator >(const Matrix& lhs, const vector rhs);
Matrix operator <=(const Matrix& lhs, const vector rhs);
Matrix operator >=(const Matrix& lhs, const vector rhs);



#endif // !DATA_CONTAINERS_H