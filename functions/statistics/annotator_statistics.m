function [sens, spec, PPV, NPV, kappa] = annotator_statistics(annotations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SENS, SPEC, PPV, NPV, KAPPA] = ANNOTATOR_STATISTICS(ANNOTATIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a set of annotations and calculates the sensitivity, 
% specificity, positive predictive value, negative predictive value, and 
% kappa relative to 50% agreement between the annotators.
%
% See our paper for details regarding this process.
% 
% Input:
%
%      annotations - a 3 dimensional matrix, n x m x p, where n and m are
%                    the image's size and p is the number of annotations. 
%                    Within which a one represents an object's location and 
%                    zero none.
%
% Output:
%
%      sens - a 1 x p vector containing the sensitivity calculated for each
%             annotation.
%
%      spec - a 1 x p vector containing the specificty calculated for each
%             annotation.
%
%      PPV - a 1 x p vector containing the positive predictive value 
%            calculated for each annotation.
%
%      NPV - a 1 x p vector containing the negative predictive value 
%            calculated for each annotation.
%
%      kappa - a 1 x p vector containing the kappa value  calculated for 
%              each annotation.
%
%	T. Lampert, A. Stumpf, and P. Gancarski, 'An Empirical Study into 
%       Annotator Agreement, Ground Truth Estimation, and Algorithm 
%       Evaluation', IEEE Transactions on Image Processing 25 (6): 
%       2557–2572, 2016.
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


gt              = sum(annotations, 3)/size(annotations,3);
m_gt            = zeros(size(gt));
m_gt(gt >= 0.5) = 1;


tp = zeros(size(annotations,3));
fp = zeros(size(annotations,3));
tn = zeros(size(annotations,3));
fn = zeros(size(annotations,3));
for i = 1:size(annotations,3)
    x = (m_gt * 2) - annotations(:,:,i);
    tp(i) = sum(x(:) == 1);
    fp(i) = sum(x(:) == -1);
    tn(i) = sum(x(:) == 0);
    fn(i) = sum(x(:) == 2);
end


sens  = zeros(size(annotations,3), 1);
spec  = zeros(size(annotations,3), 1);
PPV   = zeros(size(annotations,3), 1);
NPV   = zeros(size(annotations,3), 1);
kappa = zeros(size(annotations,3), 1);
for i = 1:size(annotations,3)
    
    sens(i)  = tp(i) / (tp(i)+fn(i));
    spec(i)  = tn(i) / (tn(i)+fp(i));
    PPV(i)   = tp(i) / (tp(i)+fp(i));
    NPV(i)   = tn(i) / (tn(i)+fn(i));
    
    p_e = (((tp(i)+fp(i))/numel(m_gt)) * ((tp(i)+fn(i))/numel(m_gt))) + (((tn(i)+fn(i)) / numel(m_gt)) * ((tn(i)+fp(i)) / numel(m_gt))); % = annotater says yes + gt says yes + annotater says no + gt says no
    kappa(i) = (((tp(i) + tn(i))/numel(m_gt)) - p_e) / (1-p_e);
    
end