// THIS SHOULD EVENTUALLY BE TURNED INTO A C++ CLASS
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>

// Main function --------------------------------------------------------------
void setCover(
    double *reconstruction,
    double *axis,
    double *radius,
    double *depth,
    double *price,
    const double *encoding,
    const double *cost,
    const double *scales,
    const double **filters,
    const double ws,
    const std::string shape,
    const size_t height,
    const size_t width,
    const size_t numChannels,
    const size_t numScales )

{
    const double BIG = 1e30;    // very big positive constant
    const size_t numImagePixels = height*width;
    size_t numPixelsCovered = 0;

    // Mark as visited pixels that cannot be covered when using disks
    if (shape == "disk") {
        for (size_t row = 0; row < height; ++row) {
            for (size_t col = 0; col < scales[0]; ++col) {
                depth[row + col*height] = depth[row + (width - col - 1)*height] = 1;
                ++numPixelsCovered;
            }
        }
        for (size_t col = 0; col < width; ++col) {
            for (size_t row = 0; row < scales[0]; ++row) {
                depth[row + col*height] = depth[height - row - 1 + col*height] = 1;
                ++numPixelsCovered;
            }
        }
    }

    // Compute how many pixels are covered by each r-disk
    std::vector<double> diskAreas(numScales);
    for (size_t scale=0; scale < numScales; ++scale) {
        for (size_t i=0; i < scale*scale; ++i) { 
            diskAreas[scale] += filters[scale][i];
        }
    }

    // Compute diskCost, costPerPixel and diskCostEffective
    std::vector<double> diskCost(cost, cost+height*width*numScales); // copy of cost
    std::vector<double> diskCostPerPixel(diskCost);
    std::vector<double> diskCostEffective(height*width*numScales);
    std::vector<double> numNewPixelsCoveredByDisk(height*width*numScales);
    for (size_t scale=0; scale < numScales; ++scale) {
        for (size_t p=0; p < numImagePixels; ++p ) {
            size_t idx = p + scale*numImagePixels;
            numNewPixelsCoveredByDisk[idx] = diskAreas[scale];
            diskCostPerPixel[idx] /= numNewPixelsCoveredByDisk[idx];
            diskCostEffective[idx] = diskCostPerPixel[idx] + ws / scales[scale];
        }
    }

    // We keep the two smallest elements in diskCostEffective and we update them
    // on the fly as we update the diskCosts, to avoid recomputing the minimum
    // at each iteration of the greedy algorithm.
    double minDiskCostEffective       = BIG;
    double secondMinDiskCostEffective = BIG;
    size_t indMinDiskCostEffective, indSecondMinDiskCostEffective;
    for (size_t i=0; i < numImagePixels*numScales; ++i) {
        if (diskCostEffective[i] < minDiskCostEffective) {
            // min1 < min2 --> newMin < min1
            secondMinDiskCostEffective = minDiskCostEffective;
            indSecondMinDiskCostEffective = indMinDiskCostEffective;
            minDiskCostEffective = diskCostEffective[i];
            indMinDiskCostEffective = i;
        }
        else if (diskCostEffective[i] < secondMinDiskCostEffective) {
            // min1 < min2 --> min1 < newMin
            secondMinDiskCostEffective = diskCostEffective[i];
            indSecondMinDiskCostEffective = i;

        }
    }

    while ( numPixelsCovered < numImagePixels) {

        // If selected cost is inf return error
        if (std::isinf(minDiskCostEffective)) {
            std::perror("Stopping: selected disk has infinite cost.")
        }

        // Get subcript indices
        const size_t row,col,scale;

        // Find NEW pixels covered by selected disk (D). If no new pixels are 
        // covered exit with error.

        // Update MAT fields

        // Update numNewPixelsCoveredByDisk, diskCostPerPixel, diskCost, 
        // diskCostEffective
    }


}