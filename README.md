# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:
- Add an additional pre-processing step that accumulates consensus from neighboring 
    disks with similar histograms at similar radii and incorporate it in the total score 
    (perhaps in the hopes of avoiding rescoring in the main while loop?).
- Assign high cost to low-scale disks but handle them as special cases towards the end, to cover leftover pixels.
- Use max(binDistance) when comparing two histograms for maximality error (and perhaps reconstruction error?).
- Use a non-linear transformation on the chi2-distance.



### QUESTIONS
- Should I compare histogram representations for the reconstruction error as well, or just for the maximality error?
- Should I adjust the B and sigma depending on r?

### Low priority
- Add mask shapes into a separate matlab class
- Fill in missing histogram distance metrics in histogramDistance()
- maybe squeeze cases with sum() and max() in histogramDistance()?
- maybe add chi2-gaussian as sub-case of chi-quadratic in histogramDistance()?

