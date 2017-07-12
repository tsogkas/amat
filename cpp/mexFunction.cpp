// Interface with MATLAB -------------------------------------------------------
#include "matrix.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // Parse inputs
    mxArray *mxMat = prhs[0];
    mxArray *mxZeroLabNormalized = prhs[1];

    // Initialize outputs
    const mxArray *mxEncoding = mxGetField(mxMat, 0, "encoding");
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
    const size_t numScales R = dims[3];
    const double *encoding   = mxGetPr(mxEncoding); 
    const double *cost       = mxGetPr(mxGetField(mxMat, 0, "cost"))
    const double *scales     = mxGetPr(mxGetField(mxMat, 0, "scales"));
    const double *filters[numScales];
    const double ws = mxGetScalar(mxGetField(mxMat, 0, "ws"));
    const mxArray *mxFilters = mxGetField(mxMat, 0, "filters");
    for (size_t r=0; r<numScales; ++r) {
        filters[r] = mxGetPr(mxGetCell(mxFilters, r)) ;
    }
    const std::string shape  = mxArrayToString(mxGetField(mxMat, 0, "shape"));



    // Do the work
    setCover(reconstruction, axis, radius, depth, price, // outputs
            encoding, cost, scales, filters, ws, shape,
            height, width, numChannels, numScales); // rest

    // Assign outputs to plhs
    plhs[0] = mxReconstruction;
    plhs[1] = mxAxis;
    plhs[2] = mxRadius;
    plhs[3] = mxDepth;
    plhs[4] = mxPrice;

    return;
}