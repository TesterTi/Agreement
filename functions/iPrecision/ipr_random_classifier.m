function [p, r, auc] = ipr_random_classifier(pi1, pi2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the performance of a random classifier in integrated 
% Precision-Recall space for specific dataset skew range, pi1 and pi1, and 
% the area under its curve.
% 
% Detailed in Appendix C, Section 4 of the paper:
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



if ~isscalar(pi1) || ~isfloat(pi1) || pi1 > 1 || pi1 < 0
    error('Pi 1 must be a floating point scalar between 0 and 1');
end

if ~isscalar(pi2) || ~isfloat(pi2) || pi2 > 1 || pi2 < 0
    error('Pi 2 must be a floating point scalar between 0 and 1');
end




% Calculate the iP-R curve
r = 0:0.001:1;
p = ones(size(r)) * (pi1 + pi2) / 2; % Eq. 21


% Calculate the AUC
auc = (pi1 + pi2) / 2;


end