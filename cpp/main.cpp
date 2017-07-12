// THIS FUNCTION IS ONLY USED FOR DEBUGGING AND DOES NOT RETURN MEANINGFUL RESULTS
#include <iostream>

int main(int argc, char *argv[]) {
    // Use dummy arguments and sizes
    const size_t height      = 256; 
    const size_t width       = 256;
    const size_t numChannels = 3;
    const size_t numScales R = 40;

    // C++ pointers to outputs (initialize to zeros)
    double reconstruction[height*width*numChannels] = {0};
    double axis[height*width*numChannels] = {0};
    double price[height*width]  = {0};
    double radius[height*width] = {0}; 
    double depth[height*width]  = {0};

    // Necessary C++ variables to keep cpp independend of the mex interface
    const double encoding[height*width*numChannels*numScales] = {0};
    const double cost[height*width*numScales] = {0};     
    const double scales[40];
    const double *filters[numScales];
    const double ws = 1e-4;
    const std::string shape = "disk";
    for (size_t r=0; r<numScales; ++r) {
        scales[r] = r+2;
        filters[r] = new double[ scales[r]*scales[r] ] ;
    }
    
    // Do the work
    setCover(reconstruction, axis, radius, depth, price, // outputs
            encoding, cost, scales, filters, ws, shape, 
            height, width, numChannels, numScales); // rest
    
    // Delete allocated filters
    for (size_t r=0; r<numScales; ++r) { delete [] filters[r] } ; 
    