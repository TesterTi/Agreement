function [corr, p] = image_feature_agreement_correlation(img, agreement)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CORR, P] = IMAGE_FEATURE_AGREEMENT_CORRELATION(IMG, ANNOTATIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes an image and the agreement levels and calculates the 
% correlation between the two.
%
% See our paper for details regarding each of this process.
% 
% Input:
%
%      img - the, n x m, image within which annotation is performed. The 
%            image should have one channel.
%
%      agreement - a 2 dimensional matrix, n x m, where n and m are
%                    the image's size. Within which each pixel value 
%                    represents the number of annotators that marked the 
%                    pixel as containing an object, such that 
%                    agreement = sum(annotations, 3).
%
% Output:
%
%      corr - a scalar representating the correlation strength.
%
%      p - a scalar representing the statistical significance of the
%          correlation (p-value).
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


% intensity correlation coefficient
img = (double(img(:,:,:))/double(max(img(:))));

[corr, p] = corrcoef(img(:), agreement(:), 'alpha', 0.01);
corr = corr(2);
p    = p(2);