function [p, r, au] = wipr_unachievable(pi1, pi2, w_function)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the unachievable area of weighted integrated Precision-Recall 
% space for a specific dataset skew range, pi1 and pi2, the weights defined 
% by the anonymous function w_function and the area under its curve.
% 
%
% Detailed in Section 3.3.2 and Appendix C, Section 3 of the paper:
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


if ~isscalar(pi1) || ~isfloat(pi1) || pi1 > 1 || pi1 < 0
    error('Pi 1 must be a floating point scalar between 0 and 1');
end

if ~isscalar(pi2) || ~isfloat(pi2) || pi2 > 1 || pi2 < 0
    error('Pi 2 must be a floating point scalar between 0 and 1');
end

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

    pi_m = (1 / (pi2-pi1)) * integral(w_function, pi1, pi2);
    
    p = zeros(size(r));
    for r_ind = 1:numel(r)
        f = @(pi) (((w_function(pi)./pi_m) .* pi .* r(r_ind)) ./ ((pi .* r(r_ind)) + (1 - pi)));
        p(r_ind) = (1/(pi2-pi1)) * integral(f, pi1, pi2); % Eq. (19)
    end

else

    % Discrete -- faster but (slightly) less accurate
    
    pi = pi1 : (pi2 - pi1) / 1000 : pi2;
    p = zeros(size(r));
    pi_n = (w_function(pi)) ./ mean(w_function(pi));
    for r_ind = 1:numel(r)
        p(r_ind) = sum((pi_n .* pi .* r(r_ind)) ./ ((pi .* r(r_ind)) + (1 - pi)));
    end
    p = p / numel(pi); % Eq. (19)

end


% Calculate the area under the curve using trapezoidal integration of 
% the P-R curve, Eq. (20)
au = auc(p,r);

end