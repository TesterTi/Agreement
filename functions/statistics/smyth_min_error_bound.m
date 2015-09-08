function error_bound = smyth_min_error_bound(annotations)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERROR_BOUND = SMYTH_MIN_ERROR_BOUND(ANNOTATIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function takes a set of annotations and calculates the theoretical 
% minimum bound on the error that can have occurred during the annotation 
% process. This is calculated according to:
%
%     P. Smyth, 'Bounds on the mean classification error rate of multiple 
%     experts' Pattern Recognition Letters 17(2): 1253-1257, 1996.
%
% Input:
%
%      annotations - a 3 dimensional matrix, n x m x p, where n and m are
%                    the image's size and p is the number of annotations. 
%                    Within which a one represents an object's location and 
%                    zero none.
%
% Output:
%
%      error_bound - a scalar representing the minimum error bound.
%
%
%	T. Lampert, A. Stumpf, and P. Gancarski, 'An Empirical Study into 
%       Annotator Agreement, Ground Truth Estimation, and Algorithm 
%       Evaluation', (submitted).
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


annotations(annotations > 0) = 1;

agreement = sum(annotations, 3);

agreement = double(agreement);

K = max(max(agreement));
N = numel(agreement);
%N = sum(gt(gt > 0));

error_bound = 0;
for i = 1:numel(agreement)
    
    error_bound = error_bound + (K - max(agreement(i), K-agreement(i)));
    
end

error_bound = error_bound / (K * N);



% Binary Case

    
K = max(max(agreement));
N = numel(agreement);

error_bound = 0;
for d = 1:K-1
    
    %n_d = sum(gt(gt == d));
    n_d = sum(sum(agreement == d));
    
    error_bound = error_bound + (n_d * min([K-d, d]));
    
end


error_bound = error_bound / (K * N);

