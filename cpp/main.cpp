// THIS FUNCTION IS ONLY USED FOR DEBUGGING AND DOES NOT RETURN MEANINGFUL RESULTS
#include <iostream>
#include <random>
#include "AMAT.h"

int main() 
{
    // Use dummy arguments and sizes
    const size_t numRows     = 256; 
    const size_t numCols     = 256;
    const size_t numChannels = 3;
    const size_t numScales   = 40;

    // C++ pointers to outputs (initialize to zeros)
    double *reconstruction = new double[numRows*numCols*numChannels]();
    double *axis    = new double[numRows*numCols*numChannels]();
    double *price   = new double[numRows*numCols]();
    double *radius  = new double[numRows*numCols](); 
    double *depth   = new double[numRows*numCols]();

    // Necessary C++ variables to keep cpp independend of the mex interface
    double *encoding = new double[numRows*numCols*numChannels*numScales]();
    double *cost = new double[numRows*numCols*numScales]();     
    const double ws = 1e-4;
    const std::string shape = "disk";
    // non const versions of radiusAtScale and filters 
    double radiusAtScale[numScales]; 
    const double *filters[numScales];

    // Initialize filters
    double *dummyFilter[numScales];
    for (int scale = 0; scale < numScales; ++scale) {
        int r = scale + 2;
        radiusAtScale[scale] = r;
        dummyFilter[scale] = new double[ (2*r+1)*(2*r+1) ]() ;
        for (int x = -r; x <= r; ++x) {
            for (int y = -r; y <= r; ++y) {
                dummyFilter[scale][y + r + (2 * r + 1) * (x + r)] = static_cast<double>((x*x + y*y) <= r*r);
            }
        }
        filters[scale] = dummyFilter[scale];
    }

    // Initialize cost
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> randomCost(0, 1);
    for (size_t i = 0; i < numRows * numCols * numScales; ++i) {
        cost[i] = randomCost(gen);
    }

    // Do the work
    setCover(reconstruction, axis, radius, depth, price,        // outputs
             encoding, cost, radiusAtScale, filters, ws, shape, // rest
             numRows, numCols, numChannels, numScales
    ); 
    
    // Free memory
    delete[] reconstruction;
    delete[] axis;
    delete[] price;
    delete[] radius;
    delete[] depth;
    delete[] encoding;
    delete[] cost;
    for (size_t r = 0; r < numScales; ++r) {
        delete[] filters[r];
    }
    return 0;
}