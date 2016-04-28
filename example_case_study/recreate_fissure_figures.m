function recreate_fissure_figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECREATE_FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function recreates the figures from the fissure case study in the paper:
%
%	T. Lampert, A. Stumpf, and P. Gancarski, 'An Empirical Study into 
%       Annotator Agreement, Ground Truth Estimation, and Algorithm 
%       Evaluation', IEEE Transactions on Image Processing 25 (6): 
%       2557–2572, 2016.
%
% This is a good starting point to follow and understand the included
% functions. All of the outputs will be placed within a folder called
% 'output'. This will include:
%
%       - various figures in eps and matlab format.
%
%       - agreement_statistics.txt containing the statistics calculated
%         using the annotations. (See agreement_analysis.m for more
%         information.)
%
%       - detector_ranks.txt containing the performance of each detector 
%         calculated using the different ground truths, and finally the 
%         rankings observed when using different ground truths. (See
%         detector_analysis.m for more information.
%
%       - a folder called gts which will contain the ground truths
%         calculated according to the annotations.
%
%
%   Copyright 2013
%
%
%   This is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This software is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this software. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup the paths
path(path, '..');
path(path, './data');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD THE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the image
img = imread('./data/image.tif');

% Load all of the annotations
annotations = load_all_annotations;

% Load all of the detections (you can replace these with any number of 
% detections that you may have)
detections(:,:,1) = imread('./data/detection_2DGWLC.tif');
detections(:,:,2) = imread('./data/detection_GAUSS.tif');
detections(:,:,3) = imread('./data/detection_CS.tif');
detections(:,:,4) = imread('./data/detection_TH.tif');

% Load the training mask (2DGWLC is the only supervised method so we use
% its training mask for all detector to make the comparison fair). Comment
% out these two lines if you are using unsupervised methods.
training_masks = imread('./data/training_mask_2DGWLC.png');
training_masks = training_masks(:,:,ones(1,4)); 

% Create the detector labels (optional)
detector_labels = {'2DWLC', 'GAUSS', 'CS', 'TH'};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM THE ANALYSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform the annotator analyses
agreement_analysis(img, annotations, detections, training_masks, detector_labels, 'ipr', 0.1, 0.5);

% Here we have pre-computed the ground truths to save time, to re-compute 
% them (or use this function for another dataset) uncomment lines 105 and
% 106 and comment the following 8 lines...
gts(:,:,1) = imread('./output/gts/gt_05.tif');
gts(:,:,2) = imread('./output/gts/gt_075.tif');
gts(:,:,3) = imread('./output/gts/gt_any.tif');
gts(:,:,4) = imread('./output/gts/gt_excl_05.tif');
gts(:,:,5) = imread('./output/gts/gt_exclude.tif');
gts(:,:,6) = imread('./output/gts/gt_lsml.tif');
gts(:,:,7) = imread('./output/gts/gt_staple.tif');
gt_labels = {'Agr. 50%', 'Agr. 75%', 'Any', 'Excl. 50%', 'Excl.', 'LSML', 'STAPLE'};

% Calculate the ground truths (this part will take a while)
%fprintf('Calculating ground truths (this will take some time)\n');
%[gts, gt_labels] = calculate_GTs(annotations);

% Perform the detector analyses (performance and ranking according to the 
% different ground truths)
detector_analysis(gts, annotations, detections, training_masks, gt_labels, detector_labels, 'ipr', 0.1, 0.5);