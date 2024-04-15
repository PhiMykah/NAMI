#include "ReadNPY.h"

using std::ios;

std::string toString(char* a, int size){
    std::string s = "";
    for (int i = 0; i < size; i++) {
        if (a[i] != '\0'){
            s = s + a[i];
        }
    }

    return s;
}

Matrix loadNPYFile(char file_path[]){
    std::ifstream array_file;
    fprintf(stderr, "Opening %s\n", file_path);

    // Open the file as a binary file
    array_file.open(file_path, ios::binary | ios::in);
    // Ensure that loading only continues if file is properly loaded
    if (!array_file){
        fprintf(stderr, "Unable to open file!\n");
        array_file.close();
        return Matrix();
    } 

    fprintf(stderr, "File opened!\n");

    // Seek to the header size and read HEADER_LEN provided by the data.
    array_file.seekg(8, ios::beg);
    char size_bytes[2];
    array_file.read(size_bytes, 2);

    // Convert bytes into little-endian unsigned short int
    // Should be multiple of 64 when 10 is added e.g. 118+10 % 64 = 0
    unsigned short h_size = ((unsigned short)size_bytes[1] << 8) | (unsigned char)size_bytes[0];

    // Initialize header char array and load header into array
    char * hdr = new char [h_size];
    array_file.read(hdr, h_size);

    // Store start of data stream position for populating matrix
    std::streampos data_start = array_file.tellg();

    // Convert header to string to allow for string parsing methods
    std::string header = toString(hdr, h_size);

    // Obtain header parameters to ensure proper data reading
    HeaderNPY headerParams = parseHeader(header, h_size);

    if (headerParams.fortran_order) {
        fprintf(stderr, "Fortran order arrays are currently unsupported.");
        array_file.close();
        return Matrix();
    }

    // Form Matrix from data stream
    return PopulateMatrix(array_file, headerParams, data_start);
}


HeaderNPY parseHeader(std::string header, int header_size){

    std::string data_code = header.substr(header.find("<")+1, 2);
    
    std::string fort_order = header.substr(header.find("'fortran_order': ")+17, 1);

    std::string shape = header.substr(header.find("(")+1, header.find("),")-header.find("(")-1);

    // Initialize variables
    DataType dtype;
    bool isFortranOrder = false;
    int M;
    int N;

    // Set datatype based on header code
    if (data_code == "f8") {dtype = DataType::f8;}
    else if (data_code == "f4") {dtype = DataType::f4;}
    else if (data_code == "i8") {dtype = DataType::i8;}
    else if (data_code == "i4") {dtype = DataType::i4;}
    else {dtype = DataType::f8;}

    // Set fortran boolean based on the first letter of the parameter
    if (fort_order == "T") {isFortranOrder = true;}

    // Set M and N
    std::string M_str = shape.substr(0, shape.find(","));
    std::string N_str = shape.substr(shape.find(", ")+2);

    M = stoi(M_str);
    N = stoi(N_str);
    
    // std::string dtype_output;
    // switch (dtype){
    // case DataType::f4:
    //     dtype_output = "f4";
    //     break;
    // case DataType::i8:
    //     dtype_output = "i8";
    //     break;
    // case DataType::i4:
    //     dtype_output = "i4";
    //     break;
    // case DataType::f8:    
    // default:
    //     dtype_output = "f8";
    //     break;
    // } 
    // std::cout << "Datatype: " << dtype_output << std::endl;
    // std::cout << "Fortran Order: ";
    // if (isFortranOrder) {std::cout << "True\n";} else {std::cout << "False\n";}
    // std::cout << "Shape: " << M << ", " << N << std::endl;

    return HeaderNPY(dtype, isFortranOrder, M, N);
}

Matrix PopulateMatrix(std::ifstream &array_file, HeaderNPY header, std::streampos data_start){
    int byte_length;

    switch (header.dtype)
    {
    case DataType::f8:
        byte_length = sizeof(double);
        break;
    case DataType::f4:
        byte_length = sizeof(float);
        break;
    default:
        fprintf(stderr, "Datatype not yet implemented!");
        array_file.close();
        fprintf(stderr, "File Closed!\n");
        return Matrix();
        break;
    }

    // Seek to end of the stream
    array_file.seekg(0, array_file.end);
    int end = array_file.tellg();

    // Return to start of data stream and calculate length
    array_file.seekg(data_start, array_file.beg);
    int start = array_file.tellg();
    int data_length = end - start;

    // Read the data stream length to a buffer
    // char * buffer = new char [data_length];
    // array_file.read(buffer,data_length);

    if (data_length/byte_length < (header.M_size * header.N_size)) {
        fprintf(stderr, "Mismatching header shape (%i,%i) with buffer size %i",
                header.M_size, header.N_size, data_length/byte_length);
        array_file.close();
        fprintf(stderr, "File Closed!\n");
        return Matrix();
    }

    // Initialize 2D vector
    vec2D newArray;

    // For loop iterates and collects each value from data stream
    // Converts 8 bytes and 4 bytes to float
    for (int i = 0; i < header.M_size; i++)
    {   
        vector newVector;
        for (int j = 0; j < header.N_size; j++)
        {
            char value[byte_length];
            array_file.read(value, byte_length);
            float newValue;
            for (int byte = 0; byte < byte_length; byte++)
            {
                printf("% 03i ", (unsigned int)value[byte]);
            }
            std::cout << std::endl;

            newVector.push_back(newValue);
        }
        newArray.push_back(newVector);
    }
    

    array_file.close();
    fprintf(stderr, "File Closed!\n");
    return Matrix();
}