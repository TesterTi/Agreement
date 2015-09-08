function gt = agreement_gt(annotations, level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SENS, SPEC, PPV, NPV, KAPPA] = ANNOTATOR_STATISTICS(ANNOTATIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the GT according to the agreement specified in 
% level (GT >= level) unless level = 0 in which case GT > level.
%
% See our paper for more details.
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
%      gt - a n x m matrix representing the resulting ground truth. Within 
%           which a one represents an object's location and zero none.
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

gt = zeros(size(annotations(:,:,1)));

if level > 0
    
    gt(agreement >= (size(annotations,3) * level)) = 1;
    
else
    
    gt(agreement > (size(annotations,3) * level)) = 1;
    
end