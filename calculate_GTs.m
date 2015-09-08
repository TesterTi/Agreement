function [gts, gt_labels] = calculate_GTs(annotations)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [GTS, GT_LABELS] = CALCULATE_GTS(ANNOTATIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a set of annotations and calculates several ground
% truth estimations based upon them (LSML, Any position marked by an 
% annotator, 50% agreement between annotators, 75% agrement between 
% annotators, STAPLE, exclude outlier annotators, and exclude outlier 
% annotators and taking 50% agreement between the remaining). See our paper 
% for references to each of these methods.
% 
% Input:
%
%      annotations - a 3 dimensional matrix n x m x p, where n and m are
%                    the image's size and p is the number of annotations.
%                    Within which a one represents an object's location and 
%                    zero none.
%
% Output:
%
%      gts - a 3 dimensional matrix n x m x 7 matrix, in which the third
%            dimension represents different ground truths. Within which a 
%            one represents an object's location and zero none.
%
%      gt_labels - a 1 x 7 cell, which contains the names for the methods
%                  used to calculate each ground truth.
%
%      This function also creates a subdirectory called 'gts' in which it
%      saves images files of each ground truth.
%
%
%	T. Lampert, A. Stumpf, and P. Gancarski, 'An Empirical Study into 
%       Annotator Agreement, Ground Truth Estimation, and Algorithm 
%       Evaluation’, (submitted).
%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputdir = './output/gts';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



path(path, genpath([fileparts(which('calculate_GTs')) filesep 'functions']));


outputdir = correct_path(outputdir);
if ~isempty(outputdir) && ~exist(outputdir, 'dir')
    mkdir(outputdir);
end


annotations(annotations > 0) = 1;




% CALCULATE GT-ANY
gt_any = agreement_gt(annotations, 0);
imwrite(gt_any, [outputdir 'gt_any.tif'], 'TIF');


% CALCULATE GT-0.5
gt_05 = agreement_gt(annotations, 0.5);
imwrite(gt_05, [outputdir 'gt_05.tif'], 'TIF');


% CALCULATE GT-0.75
gt_075 = agreement_gt(annotations, 0.75);
imwrite(gt_075, [outputdir 'gt_075.tif'], 'TIF');


% CALCULATE GT-STAPLE
gt_staple = STAPLE(annotations);
imwrite(double(gt_staple), [outputdir 'gt_staple.tif']);


% CALCULATE GT-EXCLUDE
inliers = true(1, size(annotations, 3));
[~, outliers] = cluster_annotators(annotations);
inliers(outliers) = 0;
gt_exclude = agreement_gt(annotations(:,:,inliers), 0);
imwrite(double(gt_exclude), [outputdir 'gt_exclude.tif']);


% CALCULATE GT-EXCL-0.5
gt_excl_05 = agreement_gt(annotations(:,:,inliers), 0.5);
imwrite(gt_excl_05, [outputdir 'gt_excl_05.tif'], 'TIF');


% CALCULATE GT-LSML
gt_LSML = double(LSML(annotations));
imwrite(gt_LSML, [outputdir 'gt_lsml.tif'], 'TIF');

gts        = gt_staple;
gts(:,:,2) = gt_any;
gts(:,:,3) = gt_05;
gts(:,:,4) = gt_075;
gts(:,:,5) = gt_exclude;
gts(:,:,6) = gt_excl_05;
gts(:,:,7) = gt_LSML;

gt_labels = {'STAPLE', 'Any', 'Agr. 50%', 'Agr. 75%', 'Excl.', 'Excl. 50%', 'LSML'};

end