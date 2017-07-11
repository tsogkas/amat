#include <iostream>
#include <math.h>
#include "matrix.h"
#include "mex.h"

// Main function --------------------------------------------------------------
void mexcpp(
    const mxArray *mat,
    double *reconstruction,
    double *axis,
    double *radius,
    double *depth,
    double *price)
{
    // Get some useful quantities
    const mwSize *dims = mxGetDimensions(mxGetField(mat, 0, "encoding"));
    const mwSize H = dims[0], W = dims[1], C = dims[2], R = dims[3];
    const double BIG = 1e30;
    
    // Mark as visited the pixels that cannot be covered when using disks
    std::string shape = mxArrayToString(mxGetField(mat,0,"shape"))
    if (shape == "disk") {
        double *scales = mxGetField(mat, 0, "scales");
        double r = scales[0];
        for (size_t row = 0; row < H; ++row) {
            for (size_t col = 0; col < r; ++col) {
                depth[row][col] = depth[row][W - col - 1] = 1;
            }
        }
        for (size_t col = 0; col < W; ++col) {
            for (size_t row = 0; row < r; ++row) {
                depth[row][col] = depth[H-row-1][col] = 1;
            }
        }
        
    }

    // Compute how many pixels are covered by each r-disk

    for (int r = 0; r < R; ++r)
    {
        double *filter = mxGetPr(mxGetField(mat, r, "filters"));
        double scale = mxGetScalar(mxGetField(mat, r, "scales"))
            // Min and max center coordinates of disks whose scores must be updated
            mwSize x1 = std : max(xmin - 1 - scale, 0);
        ·
            mwSize y1 = std : max(ymin - 1 - scale, 0);
        mwSize x2 = std : min(xmax - 1 + scale, W - 1);
        mwSize y2 = std : min(ymax - 1 + scale, H - 1);
        for (int x = x1; x <= x2; ++x)
        {
            · for (int y = y1; y <= y2; ++y)
            {
                // Compute 3D linear index
                mwSize idx = r * H * W + x * H + y;
                // Fidx how many newPixelsCovered are covered by current disk
                double numPixelsSubtracted = 0;
                mwSize xsoff = std : max(0, x1 - x + r);
                mwSize ysoff = std : max(0, y1 - y + r);
                mwSize xeoff = std : min(x + r - x2, x2);
                mwSize yeoff = std : min(y + r - y2, y2);
                for (int xx = x - r + xsoff; xx <= x + r + xeoff; ++xx)
                {
                    for (int yy = y - r + ysoff; yy <= y + r + yeoff; ++y)
                    {
                        if ()
                        {
                            ++numPixelsSubtracted
                        };
                    }
                }
                // ... and subtract the count from·
                numNewPixelsCovered[idx] -= numPixelsSubtracted;
                diskCost[idx] -= numPixelsSubtracted * costPerPixel[idx];
                costPerPixel[idx] = (numNewPixelsCovered[idx] == 0) ? BIG : · diskCost[idx] / numNewPixelsCovered[idx];
                diskCostEffective[idx] = costPerPixel[idx] + · mxGetScalar(maxGetField(mat, 0, "ws")) / scale;
            }
        }
    }
}

// Interface function ---------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // Parse inputs
    mxArray *mxMat = prhs[0];
    mxArray *mxZeroLabNormalized = prhs[1];

    // Create outputs
    const mwSize *dims = mxGetDimensions(mxZeroLabNormalized); // [H,W,C]
    mxArray *mxReconstruction = mxDuplicateArray(mxZeroLabNormalized);
    mxArray *mxAxis   = mxDuplicateArray(mxZeroLabNormalized);
    mxArray *mxRadius = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    mxArray *mxDepth  = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    mxArray *mxPrice  = mxCreateDoubleMatrix(dims[0], dims[1], mxREAL);
    
    // Get real pointers to outputs
    double *reconstruction = mxGetPr(mxReconstruction);
    double *axis    = mxGetPr(mxAxis);
    double *radius  = mxGetPr(mxRadius);
    double *depth   = mxGetPr(mxDepth);
    double *price   = mxGetPr(mxPrice);

    // Do the work
    mexcpp(mxMat, reconstruction, axis, radius, depth, price);

    // Assign outputs to plhs
    plhs[0] = mxReconstruction;
    plhs[1] = mxAxis;
    plhs[2] = mxRadius;
    plhs[3] = mxDepth;
    plhs[4] = mxPrice;

    return;
}