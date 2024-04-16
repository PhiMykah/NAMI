#ifndef MAIN_H
#define MAIN_H

// CPP library includes
#include <iostream>
#include <string>

// Include local header files
#include "Datatypes/DataContainers.h"
#include "Tools/BTS/BTS.h"
#include "Tools/BTS/MeanSquareDeviation.h"
#include "Tools/BTS/ExtendedComparison.h"
#include "Tools/BTS/ComplementarySimilarity.h"
#include "Tools/BTS/Medoid.h"
#include "Tools/BTS/Outlier.h"
#include "Tools/BTS/DiversitySelection.h"
#include "Tools/BTS/NewIndex.h"
#include "Tools/BTS/Align.h"
#include "FileIO/ReadNPY.h"

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int n_atoms, float (*test)(Matrix, Metric, int));

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int N, int n_atoms, float (*test)(Matrix, Metric, int, int, float, WFactor));

void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int n_atoms, Matrix (*test)(Matrix, Metric, int));

template <typename T> void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    float percent_trimmed, int n_atoms, T criterion, 
    Matrix (*test)(Matrix, float, Metric, int, T));

template <typename T> void OutputResults(
    std::string title, Matrix matrix, std::vector<Metric> metrics,
    int percentage, int n_atoms, T start,
    std::vector<int> (*test)(Matrix, int, Metric, T, int));

#endif // !MAIN_H