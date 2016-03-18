function cont = local_contrast(img, window_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONT = LOCAL_CONTRAST(IMG, WINDOW_SIZE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function takes an image and calculates its contrast within a 
% specified window size.
%
% Input:
%
%      img - a, n x m x b, image within which annotation is performed.
%            The image may be grayscale, b = 1, or colour (RGB), b = 3.
%
%      window_size - optional, a scalar, s, determining the size of the
%                    neighbourhood within which to calculate contrast, 
%                    s x s. The default value is s = 3.
%
% Output:
%
%      cont - a, n x m, image depicting the contrast.
%
%
%	T. Lampert, A. Stumpf, and P. Gancarski, 'An Empirical Study into 
%       Annotator Agreement, Ground Truth Estimation, and Algorithm 
%       Evaluation', IEEE Transactions on Image Processing (accepted).
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



if ~exist('window_size', 'var')
    window_size = 3;
end


if size(img, 3) == 3
    
    % Convert image to L*a*b* color space
    cform2lab = makecform('srgb2lab');
    LAB = applycform(img, cform2lab);
    L = LAB(:,:,1);
    
    %L = 0.2126 * img(:,:,1) + 0.7152 * img(:,:,2) + 0.0722 * img(:,:,3);
    
else
    
    % If grayscale image use that instead
    L = img;
    
end


cont = zeros(size(img(:,:,1)));
L = double(L);

for i = window_size+1:size(img,1)-window_size
    for j = window_size+1:size(img, 2)-window_size
        
        window = L(i-window_size:i+window_size , j-window_size:j+window_size);
        
        cont(i,j) = (max(max(window)) - min(min(window))) / (max(max(window)) + min(min(window)));
        
    end
end