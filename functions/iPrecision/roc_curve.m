function [tpr, fpr, au, threshold_values] = roc_curve(response, gt, interpolate)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the Recall Operating Characteristic (ROC) curve.
%
% The response should be in the format of a N x 1 or 1 x N vector,
% where N is the number of instances.
%
% The gt should be in the format of a N x 1 or 1 x N vector, the minimum 
% value is assumed to be the negative class and maximum the positive class.
%
% 
% Detailed in the paper:
%
%   T. Lampert and P. Gancarski, 'The Bane of Skew: Uncertain Ranks and 
%   Unrepresentative Precision'. In Machine Learning 97 (1-2): 5—32, 2014.
%
%
%
%   Copyright 2012
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
% Set variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
number_of_data_points = 100;       % of the ROC Curve


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if interpolate
    warning('ROC curve interpolation not yet implemented, switching off the option...');
    interpolate = 0;
end


if ~isvector(response)
    error('The algorithm''s response should be in the format of a N x 1 or 1 x N vector, where N is the number of instances.');
end

if ~isvector(gt)
    error('The gt should be in the format of a N x 1 or 1 x N vector, where N is the number of instances.');
end

if length(response) ~= length(gt)
    error('The detection and gt should be of the same length');
end

if isa(gt, 'uint8') || isa(gt, 'uint16') || isa(gt, 'uint32') || isa(gt, 'uint64')
    gt = int8(gt);
end


pos_label = int8(max(gt(:)));
neg_label = int8(min(gt(:)));

if pos_label == neg_label
    error('The ground truth must contain at least one negative and/or positive example!');
end



% Run through threshold values to evaluate the ROC curve

threshold_values = linspace(double(min(min(response))), double(max(max(response))), number_of_data_points);

tpr  = zeros(1, number_of_data_points);
fpr  = zeros(1, number_of_data_points);
tp = zeros(1, number_of_data_points);
fp = zeros(1, number_of_data_points);
tn = zeros(1, number_of_data_points);
fn = zeros(1, number_of_data_points);

for i = 1:numel(threshold_values)
    
    detection = ones(size(response), 'int8').*neg_label;
    
    detection(response >= threshold_values(i)) = pos_label;
    
    [tpr(i), fpr(i), ~, ~, tp(i), fp(i), tn(i), fn(i)] = pr_point(detection, gt, [neg_label pos_label]);
    
end



% Interpolate along the curve if desired
if interpolate
    
    tprs_n = [];
    fprs_n = [];
    
    for j = 1:number_of_data_points-1
        
        % NOT IMPLEMENTED
        [tprs1, fprs1] = roc_interpolate(tp(j+1), tp(j), fp(j+1), fp(j), fn(j));

        tprs_n = [tprs_n, p(j), tprs1(end:-1:1)];
        fprs_n = [fprs_n, r(j), fprs1(end:-1:1)];

    end
    
    tpr = tprs_n;
    fpr = fprs_n;
    
end



% Calculate the curve's AUC
au = auc(tpr(~isnan(tpr)), fpr(~isnan(fpr)));

end