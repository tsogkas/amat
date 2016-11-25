# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:

**TOP PRIORITY**
- finish coding the smirnov transform in *imageError()*. Obtain reconstruction error by creating random sample and comparing with original image patch.
- Fix histogram encoding in *patchEncoding()* and adjust code in *drawDiskOnFigureInteractive()*.
- Add support for  _histinv()_ function.
- Investigate separability of the a,b color channels in the Lab space. Can we go from a+b to a,b?

**imageError(), patchError(), drawHistGradOnFigureInteractive(), amatHistGrad()**
- handle cases when the image boundary is crossed (maybe add infinite cost?). This is _IMPORTANT_ when the background of the input image is black, in which case it matches the default (zero)padding of conv2 and can lead to erroneous high scores for disks that extend beyond the image boundaries.

**drawHistGradOnFigureInteractive()**
- Combine reconstruction error and maximality error into the visualization

