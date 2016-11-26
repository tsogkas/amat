# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:

**imageError(), patchError(), drawHistGradOnFigureInteractive(), amatHistGrad(), diskPatch()**
- handle cases when the image boundary is crossed (maybe add infinite cost?). This is _IMPORTANT_ when the background of the input image is black, in which case it matches the default (zero)padding of conv2 and can lead to erroneous high scores for disks that extend beyond the image boundaries.
- Modify *diskPatch()/patchEncoding()* so that they take into account out of border pixels?
- Consider grouping _triplets_ of Lab values instead of each channel individually.

