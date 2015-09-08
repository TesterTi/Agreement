function [data, gt] = adjust_skew(data, gt, class_labels, p_frac)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function adjusts the skew in a dataset, data, with ground truth, gt,
% to that specified by p_frac. p_frac being the proportion of positive
% instances, i.e. p_frac = N_p / N, in the resulting dataset.
%
% Class_labels is a 1 x 2 vector containing the labels of the classes, 
% i.e. [neg_label pos_label].
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



if ~isvector(data)
    error('The data should be in the format of a N x 1 or 1 x N vector, where N is the number of instances.');
end

if ~isvector(gt)
    error('The gt should be in the format of a N x 1 or 1 x N vector, where N is the number of instances.');
end

if length(data) ~= length(gt)
    error('The detection and gt should be of the same length');
end


if ~isscalar(p_frac) || ~isfloat(p_frac) || p_frac > 1 || p_frac < 0
    error('Positive fraction must be a floating point scalar between 0 and 1');
end


pos_label = class_labels(2);
neg_label = class_labels(1);


positive_count = sum(gt(:) == pos_label);


% Extract the different classes from the dataset
data_pos = data(gt == pos_label);
gt_pos = gt(gt == pos_label);
data_neg = data(gt == neg_label);
gt_neg = gt(gt == neg_label);
        

% Check whether to throw away positive of negative instances
if positive_count/numel(gt) > p_frac
    
    adjusted_positive_count = round((numel(gt) * p_frac));
    
    
    ind = randperm(numel(data_pos));
    data_pos = data_pos(ind(1:adjusted_positive_count));
    gt_pos = gt_pos(ind(1:adjusted_positive_count));
    
else
    
    if positive_count/numel(gt) < p_frac
    
        adjusted_negative_count = round((positive_count * (1/p_frac)) - positive_count);
        
        
        ind = randperm(numel(data_neg));
        data_neg = data_neg(ind(1:adjusted_negative_count));
        gt_neg = gt_neg(ind(1:adjusted_negative_count));

    end
    
end

gt = [gt_pos; gt_neg];
data = [data_pos; data_neg];


end