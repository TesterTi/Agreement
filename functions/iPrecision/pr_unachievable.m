function [p, r, au] = pr_unachievable(pi)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the unachievable area of Precision-Recall space for
% a specific dataset skew, pi, and the area under its curve.
% 
% Detailed in the paper:
%
%   K. Boyd, V. Santos Costa, J. Davis, C. Page, 'Unachievable Region in 
%   Precision-Recall Space and Its Effect on Empirical Evaluation'. In ICML, 
%   pp. (to appear), 2012.
%
% This toolbox accompanies the paper: 
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



if ~isscalar(pi) || ~isfloat(pi) || pi > 1 || pi < 0
    error('Pi (skew) must be a floating point scalar between 0 and 1');
end




% Calculate the minimum achievable P-R curve for skew pi
r = 0:0.0001:1;
p = (pi .* r) ./ ((pi .* r) + (1 - pi));



% Calculate the area of this curve
au = auc(p, r);


end