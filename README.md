# mRNA-End-to-end-Distance-Measurement-Super-resolution
* This repository contains codes for image processing of two ends of individual single-molecule mRNAs imaged by STED super-resolution microscope in two separate channels in *Drosophila* embryos. 
* Comparing with single cell, the background noise in embryo could be higher. An additional auto-adaptive threshold thus further removes the assumed background Gaussian noise after Otsu's thresholding.


## 1. Workflow
* Adaptive spot detection: Otsu's threshold and threshold from the fitted Gaussian distribution together with volume and maximum-intensity criteria.
* Nearest-neighbor matching: Assign each AF594 RNA spot to its nearest ATTO647N RNA spot using Euclidean nearest-neighbor search.
* Physical-distance conversion: Convert centroid-to-centroid distances from pixels to nanometers using the calibrated pixel size.
* Distance-based classification: Collect paired RNA signals within 200-nm and 300-nm spatial-association thresholds.

## 2. Running the codes
* Add input images in the same path of the code file or add the folder containing input images to the search path.
* Code is MATLAB files (MATLAB R2020a) and can be run on MATLAB software.

## 3. Citation
Please cite the following paper:
Chen, Jiayi, et al. "mRNA concentration–dependent translation enables rapid and sharp patterning in resource constraint Drosophila embryos." bioRxiv (2026): 2026-06.

## 4. Contact Information
- Jiayi Chen (jiayi.chen@pku.edu.cn)
- Feng Liu (liufeng@hebut.edu.cn)
