function [ps, rs] = tipr_interpolate(TP_a, TP_b, FP_a, FP_b, tn, fn, pi_function, T)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolates between two points in temporal integrated P-R space 
% according to temporal skew characteristics defined by the anonymous 
% function pi_function within the time interval [0 T] starting at skew pi_1.
% 
% The skew_function should be an anonymous function that takes time, t, as 
% an argument and contain the initial skew, pi_1, for example
%   
%   pi_1        = 0.0001;
%   pi_function = @(t)(floor_to_one(2.^t .* pi_1));
%
% do not forget that the output should be in the range [0 1].
%
%
% Detailed in Section B.1 of the paper:
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

 
use_integral = 0;       % calculate using integral or discrete (integral is 
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


pi = (TP_b + fn) / (TP_b + fn + FP_b + tn);

total_pos = TP_b+fn;



flip = 0;
if TP_b < TP_a
    
    flip = 1;
    
    t = TP_a;
    
    TP_a = TP_b;
    TP_b = t;
    
    clear t
    
end



tp = TP_a+1:TP_b-1;



% Eq. (7): first co-ordinate, R
rs = tp ./ total_pos;   



% Eq. (7): second co-ordinate, P

s = ((FP_b - FP_a) / (TP_b - TP_a));


if use_integral

   % Integral -- slower but (slightly) more accurate
   
   ps = zeros(1,numel(tp));
   for i = 1:numel(tp)
       q = @(t) (pi_function(t) ./ (pi_function(t) + ((1 - pi_function(t)) .* (pi / (1 - pi)) .* ((FP_a + s*i) ./ tp(i)))));   % Eq. (10)
       %q = @(t) (pi_function(t) * tp(i)) ./ ((pi_function(t) * tp(i)) + ((1 - pi_function(t)) .* phi .* (FP_a + s*i)));
       ps(i) = (1/T) * integral(q, 0, T);
   end
    
    
else

   % Discrete -- faster but (slightly) less accurate
   
   t = 0 : T / 10 : T;
   pi_f = pi_function(t);
   ps = 0;
   for x = 1:numel(pi_f)
       %ps = ps + (pi(x) * tp) ./ ((pi(x) * tp) + ((1 - pi(x)) * phi * (FP_a + s*[1:numel(tp)])));  
       %ps = ps + (pi_f(x) ./ (pi_f(x) + (1 - pi_f(x)) .* (pi / (1 - pi)) .* (1 + ((FP_a + s*[1:numel(tp)]) ./ tp))));  % Eq. (10)
       ps = ps + (pi_f(x) ./ (pi_f(x) + (1 - pi_f(x)) .* (pi / (1 - pi)) .* ((FP_a + s*[1:numel(tp)]) ./ tp)));  % Eq. (10)
   end
   ps = ps/numel(pi_f);

end


if flip
    ps = fliplr(ps);
    rs = fliplr(rs);
end


end