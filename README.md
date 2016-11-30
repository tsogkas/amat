# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:

**imageError(), patchError(), drawHistGradOnFigureInteractive(), amatHistGrad(), diskPatch()**
- handle cases when the image boundary is crossed (maybe add infinite cost?). This is _IMPORTANT_ when the background of the input image is black, in which case it matches the default (zero)padding of conv2 and can lead to erroneous high scores for disks that extend beyond the image boundaries.
- Consider grouping _triplets_ of Lab values instead of each channel individually.
- Vary B, dr and gaussian sigma depending on the radius.

### Low priority
- Add mask shapes into a separate matlab class
- Fill in missing histogram distance metrics in histogramDistance()
- maybe squeeze cases with sum() and max() in histogramDistance()?
- maybe add chi2-gaussian as sub-case of chi-quadratic in histogramDistance()?
