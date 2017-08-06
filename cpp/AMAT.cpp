// THIS SHOULD EVENTUALLY BE TURNED INTO A C++ CLASS
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cassert>

const double BIG = 1e30;    // very big positive constant
    
// Main function --------------------------------------------------------------
void setCover(
    double *reconstruction,
    double *axis,
    double *radius,
    double *depth,
    double *price,
    const double *encoding,
    const double *cost,
    const double radiusAtScale[],
    const double **filters,
    const double ws,
    const std::string shape,
    const size_t numRows,
    const size_t numCols,
    const size_t numChannels,
    const size_t numScales )

{
    const size_t numImagePixels = numRows*numCols;
    size_t numPixelsCovered = 0;

    // Mark as visited pixels that cannot be covered when using disks
    // (to keep the for loop simple, there are redundant assignments)
    std::vector<bool> covered(numImagePixels);
    if (shape == "disk") {
        for (int col = 0; col < radiusAtScale[0]; ++col) {
            covered[0 + col * numRows] = true;  // first row
            covered[0 + (numCols - col - 1) * numRows] = true;
            covered[numRows - 1 + col * numRows] = true;  // last row
            covered[numRows - 1 + (numCols - col - 1) * numRows] = true;
        }
        for (int row = 0; row < radiusAtScale[0]; ++row) {
            covered[row] = true;    // first col
            covered[numRows - row - 1] = true;
            covered[row + (numCols-1) * numRows] = true;    // last col
            covered[numRows - row - 1 + (numCols-1) * numRows] = true;
        }
    }
    // We start the counter of pixels covered at the number of border pixels
    for (size_t i = 0; i < numImagePixels; ++i) {
        numPixelsCovered += covered[i];
    }

    // Compute how many pixels are covered by each r-disk
    std::vector<double> diskAreas(numScales);
    for (size_t scale=0; scale < numScales; ++scale) {
        size_t filterSize = 2 * radiusAtScale[scale] + 1;
        for (size_t i = 0; i < filterSize*filterSize; ++i) {
            diskAreas[scale] += filters[scale][i];
        }
    }

    // Compute diskCost, costPerPixel and diskCostEffective
    std::vector<double> diskCost(cost, cost+numImagePixels*numScales); // copy of cost
    std::vector<double> diskCostPerPixel(numImagePixels*numScales);
    std::vector<double> diskCostEffective(numImagePixels*numScales);
    std::vector<double> numNewPixelsCoveredByDisk(numImagePixels*numScales);
    for (size_t scale = 0; scale < numScales; ++scale) {
        for (size_t p=0; p < numImagePixels; ++p ) {
            size_t idx = p + scale*numImagePixels;
            numNewPixelsCoveredByDisk[idx] = diskAreas[scale];
            diskCostPerPixel[idx]  = diskCost[idx] / numNewPixelsCoveredByDisk[idx];
            diskCostEffective[idx] = diskCostPerPixel[idx] + ws / radiusAtScale[scale];
        }
    }

    // Flags for the NEW pixels that are covered by a selected disk.
    // We preallocate the maximum memory needed for efficiency.
    std::vector<bool>   isNewPixelCovered(numImagePixels);

    // We keep the two smallest elements in diskCostEffective and we update them
    // on the fly as we update the diskCosts, to avoid recomputing the minimum
    // at each iteration of the greedy algorithm.
    double minDiskCostEffective       = BIG;
    double secondMinDiskCostEffective = BIG;
    size_t idxMinDiskCostEffective = -1, idxSecondMinDiskCostEffective = -1;
    for (size_t i=0; i < numImagePixels*numScales; ++i) {
        if (diskCostEffective[i] < minDiskCostEffective) {
            // min1 < min2 --> newMin < min1
            secondMinDiskCostEffective = minDiskCostEffective;
            idxSecondMinDiskCostEffective = idxMinDiskCostEffective;
            minDiskCostEffective = diskCostEffective[i];
            idxMinDiskCostEffective = i;
        }
        else if ( diskCostEffective[i] < secondMinDiskCostEffective ) {
            // min1 < min2 --> min1 < newMin
            secondMinDiskCostEffective = diskCostEffective[i];
            idxSecondMinDiskCostEffective = i;

        }
    }
    assert(idxMinDiskCostEffective > 0 && idxSecondMinDiskCostEffective > 0);
    assert(idxMinDiskCostEffective != idxSecondMinDiskCostEffective);

    // Start executing the greedy algorithm
    while (numPixelsCovered < numImagePixels)
    {

        std::cout << "idxMinCost: " << idxMinDiskCostEffective << std::endl;
        // If selected cost is inf return error
        if (std::isinf(minDiskCostEffective)) {
            std::perror("Selected disk has infinite cost.");
        }

        // Get subcript indices of center and radius for selected  disk
        int scale = idxMinDiskCostEffective / numImagePixels;
        int col   = idxMinDiskCostEffective % numImagePixels / numRows;
        int row   = idxMinDiskCostEffective % numImagePixels % numRows;
        assert(scale >= 0 && scale < numScales);
        assert(col >= 0 && col < numCols);
        assert(row >= 0 && row < numRows);

        
        // Find min/max coordinates of NEW pixels covered by the selected disk
        int xmin = numCols, ymin = numRows;
        int xmax = 0, ymax = 0;
        int numNewPixels = 0;
        int r = radiusAtScale[scale];
        // Reset vector before doing anything
        std::fill(isNewPixelCovered.begin(), isNewPixelCovered.end(), false);
        for (int x = col - r; x <= col + r; ++x) {
            for (int y = row - r; y <= row + r; ++y) {
                // If this is a pixel covered by the disk, increase its depth
                size_t idx = y + x * numRows;
                size_t idxFilter = (y - row + r) + (x - col + r) * (2 * r + 1);
                if (filters[scale][idxFilter]) {
                    ++depth[idx];
                    // If it has not been covered yet, that means it's a newly 
                    // covered pixel so we set flags and check its coordinates.
                    // We also update the price for all newly covered points.
                    if (!covered[idx]) {
                        ++numNewPixels;
                        isNewPixelCovered[idx] = covered[idx] = true;
                        price[idx] = minDiskCostEffective / 
                            numNewPixelsCoveredByDisk[idxMinDiskCostEffective];
                        if (x < xmin) xmin = x; 
                        if (y < ymin) ymin = y;
                        if (x > xmax) xmax = x;
                        if (y > ymax) ymax = y;
                    }
                }
            }
        }
        numPixelsCovered += numNewPixels;

        // If no new pixels are covered exit with error.
        if (numNewPixels == 0) {
            std::perror("Selected disk covers zero (0) new pixels");
        }

        // Update remaining MAT fields
        for (int c = 0; c < numChannels; ++c) {
            size_t idx = row + col * numRows + c * numImagePixels;
            axis[idx] = encoding[idx + scale * numImagePixels * numChannels];
        }
        radius[row + col * numRows] = radiusAtScale[scale];

        // Set minDiskCostEffective before updating costs
        minDiskCostEffective = secondMinDiskCostEffective;
        idxMinDiskCostEffective = idxSecondMinDiskCostEffective;
        secondMinDiskCostEffective = BIG;
        idxSecondMinDiskCostEffective = -1;

        // TODO: make this run in parallel
        // Update costs and other quantities
        // For each radius in the range...
        for (int scale = 0; scale < numScales; ++scale) {
            r = radiusAtScale[scale]; 
            //...find range of disks that are covering at least one NEW point
            // We update scores only for the disks that are entirely
            // contained inside the image domain. Disks that are not 
            // contained in the image domain have already BIG costs
            // so they will never be selected anyway. 
            int x1 = std::max(r, xmin - r);     
            int y1 = std::max(r, ymin - r);
            int x2 = std::min(xmax + r, static_cast<int>(numCols-1-r));
            int y2 = std::min(ymax + r, static_cast<int>(numRows-1-r));
            for (int x = x1; x <= x2; ++x) {
                for (int y = y1; y <= y2; ++y) {
                    // Count how many pixels are covered by other disks...
                    int numPixelsSubtracted = 0;
                    for (int xx = x - r; xx <= x + r; ++xx) {
                        for (int yy = y - r; yy <= y + r; ++yy) {
                            numPixelsSubtracted += (
                                isNewPixelCovered[yy + xx * numRows] &&
                                filters[scale][yy-y+r + (xx-x+r)*(2*r+1)] 
                            );
                        }
                    }
                    size_t idx = y + x * numRows + scale*numImagePixels;
                    assert(idx < numImagePixels * numScales);
                    //...and subtract the respective counts from those disks
                    numNewPixelsCoveredByDisk[idx] -= numPixelsSubtracted;
                    // std::cout << "#pixels subtracted at " 
                    // << (y-y1) + (x-x1)*(y2-y1+1) 
                    // << ": " << numPixelsSubtracted << std::endl; 
                    // Update diskCost, diskCostPerPixel, and diskCostEffective
                    diskCost[idx] -= numPixelsSubtracted * diskCostPerPixel[idx];
                    diskCostPerPixel[idx] =
                        (numNewPixelsCoveredByDisk[idx] == 0) ? // avoid division with 0
                        BIG : (diskCost[idx] / numNewPixelsCoveredByDisk[idx]);
                    diskCostEffective[idx] = diskCostPerPixel[idx] + ws / r;
                    // Update pair of minimum costs
                    if (diskCostEffective[idx] < minDiskCostEffective) {
                        secondMinDiskCostEffective = minDiskCostEffective;
                        idxSecondMinDiskCostEffective = idxMinDiskCostEffective;                        
                        minDiskCostEffective = diskCostEffective[idx];
                        idxMinDiskCostEffective = idx;
                    }
                    else if (diskCostEffective[idx] < secondMinDiskCostEffective &&
                             idx != idxMinDiskCostEffective) 
                    {
                        secondMinDiskCostEffective = diskCostEffective[idx];
                        idxSecondMinDiskCostEffective = idx;
                    }
                } // end for
            } // end double for loop disks covering at least one new pixel
        }   // end for loop spanning scales

        // Make sure the same point is not selected twice by assigning BIG cost.
        for (size_t scale = 0; scale < numScales; ++scale) {
            size_t idx = row + col * numRows + scale * numImagePixels;
            diskCost[idx] = diskCostEffective[idx] = BIG;
        }

    } // end while ( numPixelsCovered < numImagePixels )


}
