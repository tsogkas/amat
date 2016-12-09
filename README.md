# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:
- IMPORTANT: computation in imageError() is wrong, IF we don't explicitly 
  assign inf score to the disks that cross image boundaries. The source of the 
  error is that we use a fixed A area (area of the filter) in the simplified
  formula, even for disks that cross image boundaries. Replace with an additional
  convolution to compute the areas, or explicitly set those pixels to inf inside
  the function.
- Add an additional pre-processing step that accumulates consensus from neighboring 
    disks with similar histograms at similar radii and incorporate it in the total score 
    (perhaps in the hopes of avoiding rescoring in the main while loop?).
    It looks like this is necessary: all fully-contained disks inside a 
    MAXIMAL disk, should have very similar encodings. This maximum distance 
    of some contained disk's encoding from the maximal disk's encoding could
    be used as an maximality or reconstruction error term.
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

