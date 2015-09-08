function [p, r, auc] = wipr_random_classifier(pi1, pi2, w_function)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the performance of a random classifier in weighted integrated 
% Precision-Recall space for specific dataset skew range, pi1 and pi1, the 
% weights defined by the anonymous function w_function and the area under 
% its curve.
% 
%
% Detailed in Section 3.3.2 and Appendix C, Section 4 of the paper:
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



% Calculate the iP-R curve
r = 0:0.001:1;

if use_integral

    % Integral -- slower but (slightly) more accurate

    pi_m = (1 / (pi2-pi1)) * integral(w_function, pi1, pi2);
    q = @(pi)((w_function(pi)./pi_m) .* pi);
    p = ones(size(r)) .* ((1/(pi2-pi1))*integral(q, pi1, pi2));

else

    % Discrete -- faster but (slightly) less accurate

    pi = pi1 : (pi2 - pi1) / 1000 : pi2;
    pi_n = (w_function(pi)) ./ mean(w_function(pi));
    
    p = pi_n .* pi;
    p = sum(p(~isnan(p))) / numel(p(~isnan(p)));       % Eq. (7)

end
    
    




% Calculate the AUC
auc = (pi1 + pi2) / 2;


end