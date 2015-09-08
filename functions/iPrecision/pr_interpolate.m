function [p_i, r_i] = pr_interpolate(TP_a, TP_b, FP_a, FP_b, fn)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolates between two points in P-R space.
% 
% Detailed the paper:
%
%   J. Davis and M. Goadrich, 'The relationship between precision-recall 
%   and ROC curves'. In ICML, pp. 233-240, 2006.
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




total_pos = TP_b+fn;



flip = 0;
if TP_b < TP_a
    
    flip = 1;
    
    t = TP_a;
    
    TP_a = TP_b;
    TP_b = t;
    
    clear t
    
end



tps = TP_a+1:TP_b-1;

% First co-ordinate, R
r_i = tps ./ total_pos;   

% Second co-ordinate, P
s = ((FP_b - FP_a) / (TP_b - TP_a));
p_i = (tps ./ (tps + FP_a + s*(1:numel(tps)))); 



if flip
    p_i = fliplr(p_i);
    r_i = fliplr(r_i);
end

end