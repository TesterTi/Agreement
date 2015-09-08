function [p, r, au] = tipr_unachievable(pi_function, T)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the unachievable area of temporal integrated Precision-Recall 
% space according to the temporal skew characteristics defined by the 
% anonymous function pi_function within the time interval [0 T], and the 
% area under its curve.
% 
% Detailed in Appendix C, Section 3 of the paper:
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

 
use_integral = 1;       % calculate using integral or discrete (integral is 
                        % slower but more accurate)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_integral
    v = ver;
    if ~any(strcmp('Symbolic Math Toolbox', {v.Name}))
        use_integral = 0;
        ST = dbstack;
        warning(['Symbolic Math Toolbox not found, turning off numerical integral evaluation (results may not be as accurate as with this toolbox). To remove this warning turn off use_integral in ' ST(1).file]);
    end
    clear v ST
end



% Calculate the minimum achievable integrated P-R curve

r = 0:0.0001:1;

if use_integral

    % Integral -- slower but (slightly) more accurate

    for r_ind = 1:numel(r)
        q = @(t)((pi_function(t) .* r(r_ind)) ./ ((pi_function(t) .* r(r_ind)) + (1 - pi_function(t)))); % Eq. (12)
        p(r_ind) = (1/T) * integral(q, 0, T); % Eq. (12)
    end

else

    % Discrete -- faster but (slightly) less accurate
    
    t = 0 : T / 10 : T;
    p = zeros(size(r));
    for r_ind = 1:numel(r)
        p(r_ind) = sum((pi_function(t) .* r(r_ind)) ./ ((pi_function(t) .* r(r_ind)) + (1 - pi_function(t)))); % Eq. (12)
    end
    p = p ./ numel(t);       % Eq. (12)

end


% Calculate the area under the curve using trapezoidal integration of the 
% P-R curve, Eq. (10)

au = auc(p, r);

end