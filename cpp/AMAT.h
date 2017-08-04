// THIS SHOULD EVENTUALLY BE TURNED INTO A C++ CLASS
#ifndef _AMAT_H_
#define _AMAT_H_
#include <iostream>

void setCover(
    double *reconstruction,
    double *axis,
    double *radius,
    double *depth,
    double *price,
    const double *encoding,
    const double *cost,
    const double *radiusAtScale,
    const double **filters,
    const double ws,
    const std::string shape,
    const size_t numRows,
    const size_t numCols,
    const size_t numChannels,
    const size_t numScales 
);
    
#endif