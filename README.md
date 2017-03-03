# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:
### HIGH PRIORITY ---------------------------------------------------------
- better weigh L and color channels in grouMedialPoints.
- eliminate disconnected 1-radius pixels from the amat.
- refineMAT() has a bug that sometimes removes some of the branches.
- make sure the encodings are correct after refining.	
- run experiments on boundary detection.
- run experiments on object proposal.
- go through recent papers on (non-semantic) segmentation and setup code for 
	segmentation experiments.
- try replacing reconstructionError() from summing all errors to returning the maximum
	of the errors of all contained disks, by using the equivalent formula

### low priority ----------------------------------------------------------
- Speedup amat.
- figure out how to balance lambda, kappa (L0Smoothing) and ws (amat) parameters.
  Use default values to begin with.
  PERHAPS: setup code to chooce lambda and kappa based on how well they cover the 
  boundaries of the BSDS500 validation set. See what edge detection algorithm 
  the L0Smoothing authors use. (-->DEPENDS ON THE CHOICE OF THE EDGE DETECTION ALGORITHM).
  Maybe directly optimize lambda and kappa values wrt to reconstruction of the image
  by the mat.
- Consider squares instead of disks. How would that affect the result?
- Add constraint so the new disks cover more than a single pixel to speed up the greedy algo.
- The uniformity/reconstruction error we are using right now is not an additive function of the total 
	reconstruction (squared) error Sum(I-I')^2, so we cannot convincingly make the argument that 
	we are optimizing wrt the reconstruction quality. We should analytically derive a relationship 
	betweem the SE and the measure we are using right now (perhaps in the form of some upper bound?)
	to sell our story better.
- In the case where you use the sketchification with L0 gradient smoothing, find a way to invert
	the smoothed image to the original one (to compare reconstruction quality).	
- replace mat2mask with "double()" version. Return depth mask instead of binary.

### Very low priority
- Create AMAT class. 
- Consider changing AMAT to CMAT.
- Add mask shapes into a separate matlab class
- Fill in missing histogram distance metrics in histogramDistance()
- maybe squeeze cases with sum() and max() in histogramDistance()?
- maybe add chi2-gaussian as sub-case of chi-quadratic in histogramDistance()?

