function [tpr, fpr, p, r, tp, fp, tn, fn] = pr_point(detection, gt, class_labels)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates a point on the Precision-Recall curve.
%
% The detection should be in the format of a N x 1 or 1 x N vector,
% where N is the number of instances. It should contain the value 
% class_labels(1) for a negative detection and class_labels(2) for a 
% positive detection. Therefore class_labels is a 1 x 2 vector containing
% the labels of the classes, i.e. [neg_label, pos_label].
%
% The gt should be in the format of a N x 1 or 1 x N vector, containing 
% the value class_labels(1) for the negative class and class_labels(2) 
% for the positive class.
%
% 
% Detailed in the paper:
%
%   T. Lampert and P. Gancarski, 'The Bane of Skew: Uncertain Ranks and 
%   Unrepresentative Precision'. In Machine Learning 97 (1-2): 5—32, 2014.
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



if ~isvector(detection)
    error('The detection should be in the format of a N x 1 or 1 x N vector, where N is the number of instances.');
end

if ~isvector(gt)
    error('The gt should be in the format of a N x 1 or 1 x N vector, where N is the number of instances.');
end

if length(detection) ~= length(gt)
    error('The detection and gt should be of the same length');
end

pos_label           = 1;
neg_label           = 0;
detection(detection == class_labels(2)) = pos_label;
gt(gt == class_labels(2))               = pos_label;
detection(detection == class_labels(1)) = neg_label;
gt(gt == class_labels(1))               = neg_label;


if ~strcmpi(class(detection), class(gt))
    detection = int8(detection);
    gt = int8(gt);
end

if ~strcmpi(class(detection), class(gt))
    detection = int8(detection);
    gt = int8(gt);
end

if isa(detection, 'uint8') || isa(detection, 'uint16') || isa(detection, 'uint32') || isa(detection, 'uint64')
    detection = int8(detection);
end

if isa(gt, 'uint8') || isa(gt, 'uint16') || isa(gt, 'uint32') || isa(gt, 'uint64')
    gt = int8(gt);
end





% Number of True Positives, False Positives, True Negatives 
% and False Negatives
difference = gt - (detection*2);
tp = sum(difference == -pos_label);
fp = sum(difference == -(2*pos_label));
tn = sum(difference == neg_label);
fn = sum(difference == pos_label);



% True Positive Rate
if (tp+fn) ~= 0
    tpr = tp / (tp+fn);
else
    tpr = 1;
end



% False Positive Rate
if (fp+tn) ~= 0
    fpr = fp / (fp+tn);
else
    fpr = 0;
end



% Recall (equal to TPR)
r = tpr;



% Precision
if tp ~= 0
    p = tp / (tp + fp);
else
    p = 0;
end


end