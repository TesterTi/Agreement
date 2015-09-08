function [p_p] = p_skew_transform(p, pi, pi_p)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Transforms precision, p, calculated within a dataset having a skew of pi 
% to that which would have resulted if precision had been calculated on a 
% dataset having a skew of pi_p.
% 
% Detailed in Section 3.1 of the paper:
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
    error('Pi must be a floating point scalar between 0 and 1');
end

if ~isscalar(pi_p) || ~isfloat(pi_p) || pi_p > 1 || pi_p < 0
    error('Pi_p must be a floating point scalar between 0 and 1');
end

if ~isfloat(p)
    warning('MATLAB:nonFloatConvertToFloatVal', 'P has been cast to double, possible loss of precision');
    p = double(p);
end



p_p = (pi_p) ./ (pi_p + ((1-pi_p) .* (pi ./ (1-pi)) * ((1 ./ p)-1)));   % Eq. (3)


end