function [clustering, outliers] = cluster_annotators(annotations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CLUSTERING, OUTLIERS] = CLUSTER_ANNOTATORS(GTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function requires the statistical toolbox.
%
% This function takes a set of annotations marking objects within an image
% and calculates the pairwise (F1) distance between each annotations. It
% then clusters these distances and determines which are outliers
% (according to the mean plus one standard deviation). See our paper for
% more details.
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
%      clustering - a matrix defining the tree hierarchical clusters (see
%                   linkage for more information). To visualise the 
%                   hierarchy use dendrogram(clustering).
%
%      outliers - a 1 x y vector (where 0 <= y <= p) which gives the index 
%                 of the annotations that are deemed outliers.
%
%
%	T. Lampert, A. Stumpf, and P. Gancarski, 'An Empirical Study into 
%       Annotator Agreement, Ground Truth Estimation, and Algorithm 
%       Evaluation', (submitted).
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



data = zeros(size(annotations,3), numel(annotations(:,:,1)));

for i = 1:size(annotations,3)
    data(i,:) = reshape(annotations(:,:,i), 1, numel(annotations(:,:,i)));
end

Y = pdist(data, 'f_distance');
clustering = linkage(Y, 'ward');

% Calculate the outliers
Y = squareform(Y);
for i = 1:size(Y, 2)
	c = Y(:, i);
    
    if i > 1 && i < size(Y, 2)
        c = [c(1:i-1); c(i+1:end)]; % remove it's own distance
    else
        if i == 1
            c = c(2:end); % remove it's own distance
        end
        if i == size(Y, 2)
            c = c(1:end-1); % remove it's own distance
        end
    end
    
    M(i) = mean(c);
end

outliers = M > (mean(M) + std(M));
outliers = find(outliers);