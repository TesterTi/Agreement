# Annotator Analysis Matlab Toolbox

Thomas A. Lampert, ICube, University of Strasbourg

This work was carried out as part of the FOSTER project, which is funded by the French Research Agency (Contract 
ANR Cosinus, ANR-10-COSI-012-03-FOSTER, 2011—2014): http://foster.univ-nc.nc/

## Introduction

This toolbox accompanies the paper (a copy is included in the example_case_study subdirectory):

T. Lampert, A. Stumpf, and P. Gancarski, 'An Empirical Study into Annotator Agreement, Ground Truth Estimation, and Algorithm Evaluation’, IEEE Transactions on Image Processing 25 (6): 2557–2572, 2016.

It contains implementations of the functions described within the paper related to agreement analysis and the
evaluation of detectors using different ground truth estimation techniques. It may also be used to recreate the 
figures for the fissure case study to gain a better understanding of the method (see the QUICK START section).

It is assumed that you have a number of annotations related to the same image.

To use the toolbox's functions, simply add the toolbox directory to Matlab's path. Within the header of each 
function may be found a short description of its purpose and in which section of the paper its mathematical 
derivation can be found.

The toolbox is separated into three main functions:

1. The agreement_analysis function calculates the statistics outlined in our paper for the collection of 
annotations passed to it.
2. The calculate_GTs function calculates ground truths using the methods outlined below:
 * the LSML algorithm;
 * the agreement of any annotator;
 * the agreement of 50% of annotators;
 * the agreement of 75% of annotators;
 * the STAPLE algorithm;
 * by excluding outliers of the annotator clustering evaluation;
 * by excluding outliers and then calculating the 50% agreement level.
3. The detector_analysis function determines the detector's performance respective to each of the ground truths
passed to the function, it then ranks the detectors based upon these performances.

If you are using the Berkeley segmentation evaluation methodology then you will need to execute the build function
to compile the mex functions (if you don't know what this is then you probably don't need to use it).


## Requirements

Some parts of the toolbox (the STAPLE method to calculate ground truth) uses java code, so a virtual machine (Java 
1.5 or later) should be installed and the system's path setup so that it can be found. The default memory allocated
to the virtual machine is 2GB, so your computer should have more than this. If you have less memory, or need to
process large images, then change the memory_allocation variable in the file /functions/gt_estimation/STAPLE.m


## Quick Start

The folder 'example_case_study' contains data and functions to recreate the fissure case study figures from within 
our publication (also located within that folder). This would be a good start to understand how to use the toolbox.
To execute the function change Matlab's path to the 'example_case_study' folder and enter 'recreate_fissure_figures'.


## Included Software

This toolbox contains parts of the following software:

* ImageJ (public domain software) - http://rsbweb.nih.gov/ij/
* STAPLE (open source software) - http://www.crl.med.harvard.edu/software/STAPLE/
* ImageJ NRRD plugin (Lesser Gnu Public License v2) - http://teem.sourceforge.net/nrrd/
* Berkeley segmentation dataset evaluation functions - http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/segbench/


## Version History

v1.0 (4/10/2013) Initial Release
[![Analytics](https://ga-beacon.appspot.com/UA-112264682-3/welcome-page?pixel)](https://github.com/igrigorik/ga-beacon)
