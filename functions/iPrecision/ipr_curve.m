function [p, r, au, threshold_values] = ipr_curve(response, gt, interpolate, pi1, pi2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the integrated Precision-Recall curve according to the skews 
% pi1 and pi2.
%
% The response should be in the format of a N x 1 or 1 x N vector,
% where N is the number of instances.
%
% The gt should be in the format of a N x 1 or 1 x N vector, the minimum 
% value is assumed to be the negative class and maximum the positive class.
%
% 
% Detailed in Section 3.3.2 of the paper:
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
number_of_data_points = 100;       % of the PR Curve


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~isscalar(pi1) || ~isfloat(pi1) || pi1 > 1 || pi1 < 0
    error('Pi 1 must be a floating point scalar between 0 and 1');
end

if ~isscalar(pi2) || ~isfloat(pi2) || pi2 > 1 || pi2 < 0
    error('Pi 2 must be a floating point scalar between 0 and 1');
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




% Run through threshold values to evaluate the iP-R curve

threshold_values = linspace(double(min(min(response))), double(max(max(response))), number_of_data_points);

p  = zeros(1, number_of_data_points);
r  = zeros(1, number_of_data_points);
tp = zeros(1, number_of_data_points);
fp = zeros(1, number_of_data_points);
tn = zeros(1, number_of_data_points);
fn = zeros(1, number_of_data_points);

for i = 1:numel(threshold_values)
    
    detection = ones(size(response), 'int8').*neg_label;
    
    detection(response >= threshold_values(i)) = pos_label;
    
    [~, ~, p(i), r(i), tp(i), fp(i), tn(i), fn(i)] = ipr_point(detection, gt, [neg_label pos_label], pi1, pi2);

end

% Interpolate along the curve if desired
if interpolate
    
    ps_n = [];
    rs_n = [];
    
    for j = 1:number_of_data_points-1
        
        [ps1, rs1] = ipr_interpolate(tp(j+1), tp(j), fp(j+1), fp(j), tn(j), fn(j), pi1, pi2);

        ps_n = [ps_n, p(j), ps1(end:-1:1)];
        rs_n = [rs_n, r(j), rs1(end:-1:1)];

    end
    
    p = ps_n;
    r = rs_n;
    
end



% Calculate the curve's AUC
au = auc(p(~isnan(p)), r(~isnan(p)));

end