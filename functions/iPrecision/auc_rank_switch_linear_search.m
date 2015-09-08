function p_s = auc_rank_switch_linear_search(response_1, response_2, gt, class_labels)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Performs an exhaustive (linear) search for the points at which the rank
% of two algorithms responses response_1 and response_2, inverts with skew. 
% This function returns an empty list if a change isn't found, otherwise,
% a list of the skews at which the inversion occurrs is returned.
%
% The maximum number of inversions that CAN occur is equal to the number 
% of times the P-R curves cross (although this is a maximum and the number 
% of inversions may be less).
%
% Class_labels is a 1 x 2 vector containing the labels of the classes, 
% i.e. [neg_label pos_label].
% 
% N.B. the speed of this function would receive a big improvement if 
%      translated into C.
% 
% Detailed in Section 3.1.1 of the paper:
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


epsilon = 0.0001;       % a sufficiently small number such that:
                        %       1 - epsilon ~= 1

x = 4;                  % number of decimal places to which to find the 
                        % point of change

output = 0;             % scrolls through the current search skew on the 
                        % command line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pos_label = int8(class_labels(2));


pi = sum(sum(gt == pos_label)) / numel(gt);



% Calculate the P-R curves at the original skew
[p_1, r_1] = pr_curve(response_1(:), gt(:), 1);
[p_2, r_2] = pr_curve(response_2(:), gt(:), 1);



% Precalculate the lowest skew AUCs to allow for a more efficient search
p_1l = p_skew_transform(p_1, pi, epsilon);
p_2l = p_skew_transform(p_2, pi, epsilon);
auc_1l = auc(p_1l, r_1);
auc_2l = auc(p_2l, r_2);



% vector to store skews in
p_s = [];



% Start at a skew of epsilon iterate up to 1-epsilon, changing skew by
% the desired precision
count = 1;
if output; fprintf('Testing pi'' = '); erase_string = repmat('\b', 1, x+2); end
for pi_p = epsilon:1/(10^x):1-epsilon
    
    if output; fprintf(['%.' num2str(x) 'f'], pi_p); end
    
    % Calculate the AUCs resulting from the next skew
    p_1h = p_skew_transform(p_1, pi, pi_p);
    p_2h = p_skew_transform(p_2, pi, pi_p);
    auc_1h = auc(p_1h, r_1);
    auc_2h = auc(p_2h, r_2);
    
    
    
    % If there is a change in the sign of subtraction then a change has
    % occurred
    if sign(auc_1l - auc_2l) ~= sign(auc_1h - auc_2h)
        % Change found
        
        p_s(count) = pi_p;
        
        count = count + 1;
    end
    
    
    
    auc_1l = auc_1h;
    auc_2l = auc_2h;
    
    if output; fprintf(erase_string); end
end
if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); end

end