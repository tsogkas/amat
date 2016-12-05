# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:

### QUESTIONS
- Should I compare histogram representations for the reconstruction error as well, or just for the maximality error?
- Should I adjust the B and sigma depending on r?
- Should I use dssim on RGB or compare on the Lab space for the reconstruction error?

### Low priority
- Add mask shapes into a separate matlab class
- Fill in missing histogram distance metrics in histogramDistance()
- maybe squeeze cases with sum() and max() in histogramDistance()?
- maybe add chi2-gaussian as sub-case of chi-quadratic in histogramDistance()?

