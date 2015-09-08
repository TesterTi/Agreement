function pi_m = auc_rank_switch_test(response_1, response_2, gt, class_labels)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tests whether the rank of two algorithms, represented by their responses 
% response_1 and response_2, inverts with skew. This function returns -1 
% if no inversion occurs, otherwise the skew at which the inversion occurs 
% will be returned.
%
% This test will only work if an odd number of inversions occur within the 
% skew range (0,1), otherwise the function auc_rank_switch_linear_search
% should be used to perform an exhaustive search. The maximum number of 
% inversions that CAN occur is equal to the number of times the P-R curves
% cross (although this is a maximum and the number of inversions may be
% less).
%
% Class_labels is a 1 x 2 vector containing the labels of the classes, 
% i.e. [neg_label pos_label].
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos_label = class_labels(2);


pi = sum(sum(gt == pos_label)) / numel(gt);



% Calculate the P-R curves at the original skew
[p_1, r_1] = pr_curve(response_1(:), gt(:), 1);
[p_2, r_2] = pr_curve(response_2(:), gt(:), 1);



% Transform precisions to high positive instance skew
p_1h = p_skew_transform(p_1, pi, 1-epsilon);
p_2h = p_skew_transform(p_2, pi, 1-epsilon);
auc_1h = auc(p_1h, r_1);
auc_2h = auc(p_2h, r_2);



% Transform precisions to low positive instance skew
p_1l = p_skew_transform(p_1, pi, epsilon);
p_2l = p_skew_transform(p_2, pi, epsilon);
auc_1l = auc(p_1l, r_1);
auc_2l = auc(p_2l, r_2);



% Test for rank switch
pi_m = -1;
if sign(auc_1h - auc_2h) ~= sign(auc_1l - auc_2l)
    
    % To make things a little more efficient auc_switch_find assumes that
    % precisions have been transformed to the highest skew in the search
    % range
    pi_m = auc_rank_switch_search(p_1h, r_1, p_2h, r_2, epsilon, 1-epsilon);
    
end


end