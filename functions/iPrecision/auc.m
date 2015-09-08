function au = auc(p, r)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the area under a curve (AUC) using trapezoidal integration.
% 
% Detailed in Appendix B, Section 2 of the paper:
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



if length(p) ~= length(r)
    error('P and R must have the same number of elements');
end

if ~isfloat(p)
    warning('MATLAB:nonFloatConvertToFloatVal', 'P has been cast to double, possible loss of precision');
    p = double(p);
end

if ~isfloat(r)
    warning('MATLAB:nonFloatConvertToFloatVal', 'R has been cast to double, possible loss of precision');
    r = double(r);
end



r_1 = r(2:end);
r   = r(1:end-1);
p_1 = p(2:end);
p   = p(1:end-1);

au = sum(abs(((r_1(:) - r(:)) .* ((p_1(:) + p(:)) / 2)))); % Eq. 12


end