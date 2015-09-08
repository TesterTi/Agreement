function pi_m = auc_rank_switch_search(p_1h, r_1, p_2h, r_2, pi_l, pi_h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Searches for the point of dataset skew in which the ranking of two 
% algorithms, measured as the AUCPR, inverts. To make things slightly more
% efficient the function assumes that the precisions, p_1h, p_2h, have
% already been transformed to the highest skew in the search range, i.e.
%
%    p_1h = p_skew_transform(p_1, pi, 1-epsilon);
%    p_2h = p_skew_transform(p_2, pi, 1-epsilon);
%
% where 1-epsilon is the upper bound on the search range and pi is the skew
% at which precision was calculated. For more information see the 
% auc_switch function.
%
% This function will only return one point at which the rank inverts, if 
% more exist within the skew range (0,1) the function 
% auc_rank_switch_linear_search should be used to perform an exhaustive 
% search. The maximum number of inversions that CAN occur is equal to the 
% number of times the P-R curves cross (although this is a maximum and the 
% number of inversions may be less).
%
% Detailed in Section 3.1.1, Algorithm 1, of the paper:
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


x = 10;                 % number of decimal places to which to find the 
                        % point of change

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check whether point of change has been found to x decimal places
if floor(10^x * pi_l) == floor(10^x * pi_h)
    
    pi_m = floor(10^x * pi_h) / 10^x; % round to x decimal places
    
    return
    
end



% Calculate the middle skew in the range
pi_m = (pi_l + pi_h) / 2;



% Transform precisions to the middle positive instance skew
p_1m = p_skew_transform(p_1h, pi_h, pi_m);
p_2m = p_skew_transform(p_2h, pi_h, pi_m);



% Calculate AUCs
auc_1h = auc(p_1h, r_1);
auc_2h = auc(p_2h, r_2);
auc_1m = auc(p_1m, r_1);
auc_2m = auc(p_2m, r_2);



if sign(auc_1m - auc_2m) ~= sign(auc_1h - auc_2h)
    
    % Change in upper half
    
    pi_m = auc_rank_switch_search(p_1h, r_1, p_2h, r_2, pi_m, pi_h);
    
else
    
    % If it's not in the upper half transform precision to the lower skew
    p_1l = p_skew_transform(p_1h, pi_h, pi_l);
    p_2l = p_skew_transform(p_2h, pi_h, pi_l);
    
    
    
    % Calculate AUCs
    auc_1l = auc(p_1l, r_1);
    auc_2l = auc(p_2l, r_2);
    
    
    
    if sign(auc_1l - auc_2l) ~= sign(auc_1m - auc_2m)
        
        % Change in lower half
        
        pi_m = auc_rank_switch_search(p_1m, r_1, p_2m, r_2, pi_l, pi_m);
        
    end
    
end


end