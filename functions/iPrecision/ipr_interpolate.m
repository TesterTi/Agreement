function [ps, rs] = ipr_interpolate(TP_a, TP_b, FP_a, FP_b, tn, fn, pi1, pi2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolates between two points in integrated P-R space.
% 
% Detailed in Appendix C, Section 1 of the paper:
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


phi = (TP_b+fn) / (FP_b+tn);

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
      p = @(pi) (pi * tp(i)) ./ ((pi * tp(i)) + ((1 - pi) .* phi .* (FP_a + s*i)));
      ps(i) = (1/(pi2-pi1)) * integral(p, pi1, pi2);
   end

else

   % Discrete -- faster but (slightly) less accurate
   
   pi = pi1 : (pi2 - pi1) / 1000 : pi2; % Skews to integrate
   ps = 0;
   for x = 1:numel(pi)
       ps = ps + (pi(x) * tp) ./ ((pi(x) * tp) + ((1 - pi(x)) * phi * (FP_a + s*[1:numel(tp)])));  
   end
   ps = ps/numel(pi);

end

if flip
    ps = fliplr(ps);
    rs = fliplr(rs);
end


end