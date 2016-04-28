function dist = f_distance(d1, d)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIST = F_DISTANCE(D1, D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function take a vector d1 and a matrix d and calculates the f1 
% distance between the vector and each row of the matrix.
%
% See our paper for details regarding each of this process and 'help pdist'
% for more information regarding the form of the function. 
%
% Input:
%
%      d1 - a 1 x N 'sample' vector 
%           
%      d - a M x N 'target' matrix
%
% Output:
%
%      dist - a M x 1 vector representing the distances between d1 and each
%             row of the matrix d.
%
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


dist = zeros(1, size(d, 1));
for i = 1:size(d,1)
    
    d2 = d(i,:);
    
    d1(d1 > 0) = 1;
    d2(d2 > 0) = 1;

    temp = d1 - 2*d2;

    TP = sum(temp == -1);
    FP = sum(temp ==  1);
    P  = sum(d2);
    
    if TP+FP > 0
        Pr = TP / (TP + FP);
    else
        Pr = 1;
    end
    
    if P > 0
        Re = TP / P;
    else
        Re = 1;
    end
    
    if Pr + Re > 0
        dist(i) = 1 - (2 * ((Pr*Re) / (Pr+Re)));
    else
        if P > 0
            dist(i) = 1 - (2 * ((Pr*Re) / (Pr+Re)));
        else
            dist(i) = 0;
        end
    end
    
end