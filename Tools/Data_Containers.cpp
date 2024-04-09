#include "Data_Containers.h"

/*
A basic quicksort function implementation

Parameters
----------
vec : vector
    Vector to sort
l_index : int
    Left-most index
r_index : int
    Right-most index
*/
void quicksort(vector &vec, int l_index, int r_index){
    int i, j, m_index, pivot;
    i = l_index;
    j = r_index;
    m_index = l_index + (r_index - l_index) / 2;
    pivot = vec[m_index];

    while(i < r_index || j > l_index){
        while (vec[i] < pivot) {
            i++;
        }
        while (vec[j] > pivot) {
            j--;
        }

        if (i <= j) {
            std::swap(vec[i],vec[j]);
            i++;
            j--;

        } else {
            if (i < r_index) {
                quicksort(vec, i, r_index);
            }
            if (j > l_index) {
                quicksort(vec, l_index, j);
            }

            return;
        }
    }
}


/*
A basic quicksort function implementation

Parameters
----------
vec : vector
    Vector to sort
indx_vector : std::vector<int>
    Vector of vec's indicies to track sorting changes
l_index : int
    Left-most index
r_index : int
    Right-most index
*/
void quicksort(vector &vec, std::vector<int> &indx_vector, int l_index, int r_index){
    int i, j, m_index;
    float pivot;
    i = l_index;
    j = r_index;
    m_index = l_index + (r_index - l_index) / 2;
    pivot = vec[m_index];

    while(i < r_index || j > l_index){
        while (vec[i] < pivot) {
            i++;
        }
        while (vec[j] > pivot) {
            j--;
        }

        if (i <= j) {
            std::swap(vec[i],vec[j]);
            std::swap(indx_vector[i],indx_vector[j]);
            i++;
            j--;

        } else {
            if (i < r_index) {
                quicksort(vec, indx_vector, i, r_index);
            }
            if (j > l_index) {
                quicksort(vec, indx_vector, l_index, j);
            }

            return;
        }
    }
}

Matrix::Matrix()
{
    this->m_array = vec2D {}; 
    this->N = 0;
    this->M = 0; 
    //this->Flatten();
}

Matrix::Matrix(vec2D array){
    this->SetArray(array);
    this->Flatten();
}

Matrix::Matrix(vec2D array, vector flat){
    this->m_array = array;
    this->N = array.size();
    if (this->N >= 1) {
        this->M = array.at(0).size();
    }
    else{
        this->M = 0;
    }
    this->m_flat = flat;
}

void Matrix::SetArray(vec2D array){
    this->m_array = array; 
    this->N = array.size();
    if (this->N >= 1) {
        this->M = array.at(0).size();
    }
    else{
        this->M = 0;
    }
    this->Flatten();
}

vec2D Matrix::GetArray(){
    return this->m_array;
}

// This doesn't check for invalid index so use wisely
void Matrix::erase(int index){
    this->m_array.erase(m_array.begin()+index);

    this->N = m_array.size();
    this->Flatten();
}

// This doesn't check for invalid index so use wisely
void Matrix::erase(std::vector<int> indicies){
    // Remove the values at those indices of the original matrix
    for (int i = 0; i < indicies.size(); i++){
        this->m_array.erase(m_array.begin()+indicies[i]);
    }

    this->N = m_array.size();
    this->Flatten();
}

void Matrix::Flatten(){
    vector flatVec(this->m_array.size()*this->m_array[0].size());
    for (int row = 0; row < this->N; row++){
        for (int col = 0; col < this->M; col++){
            flatVec[col + row * this->m_array.size()] = this->m_array[row][col];
        }        
    }
    this->m_flat = flatVec;
}

// Perform an elementwise operation on the data
Matrix Matrix::Elementwise(float (*func)(float)) const{
    vec2D newArray;
    
    for (int row = 0; row < this->N; row ++){   
        // Initialize new vector
        vector vec;

        // Square the values of the column
        for (int col = 0; col < this->M; col++){
            vec.push_back(func(this->m_array[row][col]));
        }
        // Add the new float vector to the vector container
        newArray.push_back(vec);
    }
    
    Matrix mat(newArray);

    return mat;
}

// Perform an elementwise operation on the data using target function and input value
Matrix Matrix::Elementwise(Matrix x, float (*func)(float, float)) const{
    vec2D newArray;
    
    for (int row = 0; row < this->N; row ++){   
        // Initialize new vector
        vector vec;

        // Square the values of the column
        for (int col = 0; col < this->M; col++){
            vec.push_back(func(this->m_array[row][col], x[row][col]));
        }
        // Add the new float vector to the vector container
        newArray.push_back(vec);
    }
    
    Matrix mat(newArray);

    return mat;
}

// Summate all of the matrix's values overall
float Matrix::Sum(){
    float sum = 0;
    for (int i = 0; i < this->m_flat.size(); i++)
    {
        sum += this->m_flat[i];
    }
    return sum;
}

vector Matrix::Sum(AXIS axis){
    // How far apart each value is in the flattened array
    int step;
    // How far each vector start is in the flattened array
    int increment;
    // How many points are in each vector
    int num_points;
    // How many vectors are summed
    int num_vectors;

    // Column length
    int M = this->M;

    // Row length
    int N = this->N;
    
    // Change above parameters based on which axis is chosen
    switch (axis){
    
    case AXIS::ROW:
        step = 1;
        increment = M;
        num_points = M;
        num_vectors = N;
        break;

    default:
    case AXIS::COLUMN:
        step = M;
        increment = 1;
        num_points = N;
        num_vectors = M;
        break;
    }

    // Float vector containing all of the summations
    vector Sums(num_vectors);

    // Loop over each element in the vector, sum the elements, then move to next vector
    for (int vidx = 0; vidx < increment * num_vectors; vidx += increment) {
        int sum = 0;
        for (int idx = 0; idx < num_points; idx++) {
            sum += this->m_flat[vidx + idx * step];
        }
        Sums[vidx/increment] = sum;
    }

    return Sums;

}

Matrix Matrix::Absolute(){
    return this->Elementwise([](float x){if (x < 0) {return -1*x;} else {return x;}});
}

Matrix Matrix::pow(int power){
    return this->Elementwise(cast(power, *this), POW);
}

// Add another vector to the matrix
void Matrix::push_back(vector x) {
    if (this->M == 0) {
        this->M = x.size();
    }
    else if (x.size() != this->M) {
        fprintf(stderr, "Internal Error: Vector of legnth %li cannot be added to matrix of dimensions %i, %i.\n",
                x.size(), this->N, this->M);
        return;
    }

    this->m_array.push_back(x);
    this->N += 1;
    
    for (int i = 0; i < this->M; i++){
        this->m_flat.push_back(x[i]);
    }
}

void Matrix::print(){
    for (int row = 0; row < this->N; row ++)
    {
        for (int col = 0; col < this->M; col++)
        {
            printf("% 7.2f ", this->m_array[row][col]);
        }
        std::cout << std::endl;
    }
}

void Matrix::clear(){
    this->m_array = vec2D {vector {}};
    this->N = 0;
    this->M = 0;
    this->m_flat = vector {};
}

Matrix::~Matrix(){
    this->clear();
}

// *************************
// * Reshaping and Casting *
// *************************

Matrix::operator vector() const{
    if (this->N == 0){
        return vector {};
    }

    return this->m_flat;
}

// Cast Constant value to Matrix
Matrix cast(float x, const Matrix& source_matrix){
    int N = source_matrix.N;
    int M = source_matrix.M;

    // Add N vectors
    vec2D newArray;
    for (int row = 0; row < N; row++) {
        vector newVector;
        for (int col = 0; col < M; col++) {
            newVector.push_back(x);
        }
        newArray.push_back(newVector);
    }

    Matrix mat(newArray);

    return mat; 
}


// Cast vector of floats to a 2D vector
Matrix cast(vector vec, const Matrix& source_matrix){
    int N = source_matrix.N;
    int M = source_matrix.M;

    // Fill vector with more data if less than desired N*M size
    if (vec.size() < M) {
        int v_size = vec.size();
        int i = 0;
        while (vec.size() < M) {
            vec.push_back(vec[i]);
            i = (i + 1) % v_size;
        }
    } else if (vec.size() > M) {
        vector newVec;
        for (int i = 0; i < M; i++) {
            newVec.push_back(vec[i]);
        }
        vec = newVec;
    }
    
    // Add N vectors
    vec2D newArray;
    for (int row = 0; row < N; row++) {
        newArray.push_back(vec);
    }

    Matrix mat(newArray);

    return mat;  
}


// Reshape Vector to Matrix
Matrix vectorReshape(vector vec, const Matrix& source_matrix){
    int N = source_matrix.N;
    int M = source_matrix.M;

    // Reshape from N*M vector to N x M matrix
    if (vec.size() == N*M) {
        vec2D newArray;
        for (int row = 0; row < N; row++) {

            vector newVector;
            for (int col = 0; col < M; col++){

                newVector.push_back(vec[row * N + col]);
            }
            newArray.push_back(newVector);            
        }
        Matrix mat(newArray, vec);

        return mat;
    }

    return cast(vec, source_matrix);
}


// ********************
// Overloaded Operators
// ********************

vector const& Matrix::operator[](int idx) const {
    if (idx < 0 || idx >= this->N) {
        fprintf(stderr, "Internal Error: Position %i out of range for Argument List.\n",idx);
        return EMPTY_VECTOR;
    }
    return this->m_array.at(idx);
}

Matrix Matrix::operator[](Matrix index_matrix) const {
    if ((index_matrix.N != this->N) || (index_matrix.M != this->M)) {
        fprintf(stderr, "Internal Error: Index matrix of size (%i,%i) does not match matrix of size (%i,%i)\n",
                index_matrix.N, index_matrix.M, this->N, this->M);
        return *this;
    }
    vector vec;
    for (int row = 0; row < this->N; row++) {
        for (int col = 0; col < this->M; col++){
            // Add to vector if corresponding row and column index is has a truth value
            if ((bool) index_matrix[row][col]){
                vec.push_back(this->m_array[row][col]);
            }
        }
    }
    vec2D mat = {vec};
    Matrix newMatrix(mat);
    return newMatrix;
}

// ************
// * ADDITION *
// ************

// Adds two matrices together elementwise, assumes they have the same shape
Matrix Matrix::operator+(Matrix x){
    return this->Elementwise(x, ADD);
}

// Add constant (Matrix + constant)
Matrix Matrix::operator+(float x){
    return this->Elementwise(cast(x, *this), ADD);
}

// Add constant (constant + Matrix)
Matrix operator+(float x, const Matrix& y){
    return y.Elementwise(cast(x, y), ADD);
}

// Add Vector (Matrix + vector)
Matrix Matrix::operator+(vector x){
    return this->Elementwise(cast(x, *this), ADD);
}

// Add Vector (vector + Matrix)
Matrix operator+(vector x, const Matrix& y){
    return y.Elementwise(cast(x, y), ADD);
}

// Vector-Vector addition
vector operator+(vector x, vector y){
    vector newVector;
    // Return x in case of mismatch vectors
    if (x.size() != y.size()) {
        return x;
    }

    for (int i = 0; i < x.size(); i++) {
        newVector.push_back(x[i] + y[i]);
    }

    return newVector;
}


// **************
// * SUBTRATION *
// **************

// Subtract two matrices elementwise, assumes they have the same shape
Matrix Matrix::operator-(Matrix x){
    return this->Elementwise(x, SUB);
}

// Subtract constant (Matrix - constant)
Matrix Matrix::operator-(float x){
    return this->Elementwise(cast(x, *this), SUB);
}

// Subtract constant (constant - Matrix);
Matrix operator-(float x, const Matrix& y){
    return cast(x, y).Elementwise(y, SUB);
}

// Subtract vector (Matrix - vector)
Matrix Matrix::operator-(vector x){
    return this->Elementwise(cast(x, *this), SUB);
}

// Subtract vector (vector - Matrix)
Matrix operator-(vector x, const Matrix& y){
    return cast(x, y).Elementwise(y, SUB);
}

// Vector-vector subtraction
vector operator-(vector x, vector y){
    vector newVector;
    // Return x in case of mismatch vectors
    if (x.size() != y.size()) {
        return x;
    }

    for (int i = 0; i < x.size(); i++) {
        newVector.push_back(x[i] - y[i]);
    }

    return newVector;
}

// ******************
// * MULTIPLICATION *
// ******************

// Multiply two matrices elementwise, assumes they have the same shape
Matrix Matrix::operator*(Matrix x){
    return this->Elementwise(x, SUB);
}

// Multiply constant ( matrix x constant )
Matrix Matrix::operator* (float x){
    return this->Elementwise(cast(x, *this), MULT);
}

// Multiply constant ( matrix x constant )
Matrix operator* (float x, const Matrix& y){
    return y.Elementwise(cast(x, y), MULT);
}

// Multiply vector elementwise ( matrix x vector )
Matrix Matrix::operator* (vector x){
    return this->Elementwise(cast(x, *this), MULT);
}

// Multiply vector elementwise ( vector x matrix )
Matrix operator* (vector x, const Matrix& y){
    return y.Elementwise(cast(x, y), MULT);
}


// ************
// * DIVISION *
// ************

// Divide two matrices elementwise, assumes they have the same shape
Matrix Matrix::operator/(Matrix x){
    return this->Elementwise(x, DIV);
}

// Division of matrix by constant ( matrix / constant )
Matrix Matrix::operator/(float x){
    return this->Elementwise(cast(x, *this), DIV);
}

// Divison of constant by matrix elementwise ( constant / matrix ) 
Matrix operator/(float x, const Matrix& y){
    return y.Elementwise(cast(x, y), DIV);
}

// Division of matrix by vector ( matrix / vector )
Matrix Matrix::operator/(vector x){
    return this->Elementwise(cast(x, *this), DIV);
}

// Divison of vector by matrix elementwise ( vector / matrix ) 
Matrix operator/(vector x, const Matrix& y){
    return y.Elementwise(cast(x, y), DIV);
}

// **********
// * MODULO *
// **********

Matrix Matrix::operator%(Matrix x){
    return this->Elementwise(x, MOD);
}

// Matrix modulo constant elementwise ( matrix % constant )
Matrix Matrix::operator%(float x){
    return this->Elementwise(cast(x, *this), MOD);
}

// Matrix modulo vector elementwise ( matrix % vector )
Matrix Matrix::operator%(vector x){
    return this->Elementwise(cast(x, *this), MOD);
}

// ****************
// * INEQUALITIES *
// ****************

// All values less than float value
Matrix operator <(const Matrix& lhs, const float rhs){
    return lhs.Elementwise(cast(rhs, lhs), LESS_THAN);
}

// All values greater than float value
Matrix operator >(const Matrix& lhs, const float rhs){
    return lhs.Elementwise(cast(rhs, lhs), GREATER_THAN);
}

// All values less than or equal to float value
Matrix operator <=(const Matrix& lhs, const float rhs){
    return lhs.Elementwise(cast(rhs, lhs), LEQ);
}

// All values greater than or equal to float value
Matrix operator >=(const Matrix& lhs, const float rhs){
    return lhs.Elementwise(cast(rhs, lhs), GEQ);
}

// **********
// * Powers *
// **********

// Raise Matrix to Matrix elementwise
Matrix pow(const Matrix& x, const Matrix& y){
    return x.Elementwise(y, POW);
}

// Raise Matrix to element
Matrix pow(const Matrix& x, float y){
    return x.Elementwise(cast(y, x), POW);
}

// Raise Matrix to vector
Matrix pow(const Matrix& x, vector y){
    return x.Elementwise(cast(y, x), POW);
}

// Raise value to matrix elements
Matrix pow(float x, const Matrix& y){
    return y.Elementwise(cast(x, y), POW);
}
// Raise vector to matrix elements
Matrix pow(vector x, const Matrix& y){
    return y.Elementwise(cast(x, y), POW);
}

// Raise vector to constant
vector pow(vector x, float y){
    vector newVector;
    for (int i = 0; i < x.size(); i++) {
        newVector.push_back(std::pow(x[i],y));
    }

    return newVector;
}