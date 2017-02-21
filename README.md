# amat
Code for the Appearance-MAT project

Created by Stavros Tsogkas at the University of Toronto.

### License

This code is released under the MIT License (refer to the LICENSE file for details).

### TODO:
### HIGH PRIORITY ---------------------------------------------------------
- add computation of ods,ois etc for methods with binary outputs (single threshold).
- find a way to reconstruct images for MIL and deepskel methods.
- finish testReconstruction().
- figure out how to balance lambda, kappa (L0Smoothing) and ws (amat) parameters.

### low priority ----------------------------------------------------------
- Add an additional pre-processing step that accumulates consensus from neighboring 
    disks with similar histograms at similar radii and incorporate it in the total score 
    (perhaps in the hopes of avoiding rescoring in the main while loop?).
    It looks like this is necessary: all fully-contained disks inside a 
    MAXIMAL disk, should have very similar encodings. This maximum distance 
    of some contained disk's encoding from the maximal disk's encoding could
    be used as an maximality or reconstruction error term.
    NOTE: rescoring does not solve this problem because it reduces score for 
    all contained disks in the selected disk, so smaller disks that are close 
    to the boundary are favored. We must explicitly favor the selection of disks
    whose centers are neighbors to the currently selected center.
- Assign high cost to low-scale disks but handle them as special cases towards the end, to cover leftover pixels.
- Add constraint so the new disks cover more than a single pixel to speed up the greedy algo.
- The uniformity/reconstruction error we are using right now is not an additive function of the total 
	reconstruction (squared) error Sum(I-I')^2, so we cannot convincingly make the argument that 
	we are optimizing wrt the reconstruction quality. We should analytically derive a relationship 
	betweem the SE and the measure we are using right now (perhaps in the form of some upper bound?)
	to sell our story better.
- In the case where you use the sketchification with L0 gradient smoothing, find a way to invert
	the smoothed image to the original one (to compare reconstruction quality).	


### QUESTIONS
- Should I distribute the cost of a disk over the number of its pixels before adding the maximality cost or after?
= If the disk cost is updated when some of its pixels are covered, then the division should happen BEFORE
    adding the maximality cost. Otherwise, if we manage to assign a fixed cost to each disk once in the beginning
    and then just apply the greedy algorithm, the division could go after adding the maximality cost.    

### Experiments/Applications
- For natural images use first a method that sketchifies the image and apply the algorithm on the result. See:
	[1] Image Smoothing via L_0 Gradient Minimization, Xu,Lu,Xu,Jia
	[2] Structure Extraction from Texture via Relative Total Variation, Xu,Yan,Xia,Jia
- Then use a learning method to learn from BSDS500 extracted "groundtruth" skeletons which how to combine
	histogram features for scoring "medialness" of disks. For this part however we have to also consider 
	how to locally reconstruct the image.
- Boundary extraction from the medial axis representation. Strength of the boundary may depend on its "depth",
	where depth is the number of disks the point is covered by.
- Test performance by examining _recall_ of skeleton points extracted from segments on the BSDS500 and MSCOCO datasets.
- Foreground/background separation (one simple thing is to consider as background structures that 'touch' the image border).
- Salient region/object proposal based on symmetry responses (grouped or ungrouped).
- Image summarization and application in image retrieval, or image co-segmentation.


### Very low priority
- Add mask shapes into a separate matlab class
- Fill in missing histogram distance metrics in histogramDistance()
- maybe squeeze cases with sum() and max() in histogramDistance()?
- maybe add chi2-gaussian as sub-case of chi-quadratic in histogramDistance()?

