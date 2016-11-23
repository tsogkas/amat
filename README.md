# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:
**amatSetCoverGreedy**
- Maybe add an extra regularization term that discourages disks with radius that does not agree with the radii of neighboring/enclosed disks --> this can probably be implicitly modeled by the next point in the list

**imageError(), drawHistGradOnFigureInteractive(), amatHistGrad()**
- handle cases when the image boundary is crossed (maybe add infinite cost?). This is _IMPORTANT_ when the background of the input image is black, in which case it matches the default (zero)padding of conv2 and can lead to erroneous high scores for disks that extend beyond the image boundaries.
- Change the space on which we perform comparisons to the space of histograms (penalize dissimilar histograms instead of penalizing bad raw image reconstructions).

**imageEncoding(), patchEncoding()**
- Change encoding for histograms: choosing maximum vote for each bin _independently_ is not the right way to do it. We must look for more frequent Lab _triplets_.

**drawHistGradOnFigureInteractive()**
- Combine reconstruction error and maximality error into the visualization
