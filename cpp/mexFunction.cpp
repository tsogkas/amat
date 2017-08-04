// Interface with MATLAB -------------------------------------------------------
#include "mex.h"
#include "matrix.h"
#include "AMAT.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // Parse inputs
    const mxArray *mxMat = prhs[0];
    const mxArray *mxZeroLabNormalized = prhs[1];

    // Initialize outputs
    const mxArray *mxEncoding = mxGetProperty(mxMat, 0, "encoding");
    const mwSize *dims = mxGetDimensions(mxEncoding); // [H,W,C,R]
    mxArray *mxReconstruction = mxDuplicateArray(mxZeroLabNormalized);
    mxArray *mxAxis   = mxDuplicateArray(mxZeroLabNormalized);
    mxArray *mxRadius = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    mxArray *mxDepth  = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    mxArray *mxPrice  = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    
    // C++ pointers to outputs
    double *reconstruction = mxGetPr(mxReconstruction);
    double *axis    = mxGetPr(mxAxis);
    double *radius  = mxGetPr(mxRadius);
    double *depth   = mxGetPr(mxDepth);
    double *price   = mxGetPr(mxPrice);

    // Necessary C++ variables to keep cpp independend of the mex interface
    const size_t height      = dims[0]; 
    const size_t width       = dims[1];
    const size_t numChannels = dims[2];
    const size_t numScales   = dims[3];
    const double *encoding   = mxGetPr(mxEncoding);
    const double *cost       = mxGetPr(mxGetProperty(mxMat, 0, "cost"));
    const double *scales     = mxGetPr(mxGetProperty(mxMat, 0, "scales"));
    const double *filters[numScales];
    const double ws = mxGetScalar(mxGetProperty(mxMat, 0, "ws"));
    for (size_t r=0; r<numScales; ++r) {
        filters[r] = mxGetPr(mxGetCell(mxGetProperty(mxMat, 0, "filters"), r)) ;
    }
    const std::string shape = mxArrayToString(mxGetProperty(mxMat, 0, "shape"));



    // Do the work
    setCover(reconstruction, axis, radius, depth, price, // outputs
        encoding, cost, scales, filters, ws, shape,      // rest
        height, width, numChannels, numScales
    ); 

    // Assign outputs to plhs
    plhs[0] = mxReconstruction;
    plhs[1] = mxAxis;
    plhs[2] = mxRadius;
    plhs[3] = mxDepth;
    plhs[4] = mxPrice;

    return;
}