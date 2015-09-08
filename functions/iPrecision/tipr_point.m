function [tpr, fpr, p, r, tp, fp, tn, fn] = tipr_point(detection, gt, class_labels, pi_function, T)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates a point on the tempral integrated Precision-Recall curve 
% according to temporal skew characteristics defined by the anonymous 
% function pi_function within the time interval [0 T] starting at skew pi_1.
%
% The detection should be in the format of a N x 1 or 1 x N vector,
% where N is the number of instances. It should contain the value 
% class_labels(1) for a negative detection and class_labels(2) for a 
% positive detection. Therefore class_labels is a 1 x 2 vector containing
% the labels of the classes, i.e. [neg_label, pos_label].
%
% The gt should be in the format of a N x 1 or 1 x N vector, containing 
% the value class_labels(1) for the negative class and class_labels(2) 
% for the positive class.
%
% The skew_function should be an anonymous function that takes time, t, as 
% an argument and contain the initial skew, pi_1, for example
%   
%   pi_1        = 0.0001;
%   pi_function = @(t)(floor_to_one(2.^t .* pi_1));
%
% do not forget that the output should be in the range [0 1].
% 
% 
% Detailed in Section 3.2 of the paper:
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

 
use_integral = 0;       % calculate using integral or discrete (integral is 
                        % slower but more accurate)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~isvector(detection)
    error('The detection should be in the format of a N x 1 or 1 x N vector, where N is the number of instances.');
end

if ~isvector(gt)
    error('The gt should be in the format of a N x 1 or 1 x N vector, where N is the number of instances.');
end

if length(detection) ~= length(gt)
    error('The detection and gt should be of the same length');
end

if ~strcmpi(class(detection), class(gt))
    detection = int8(detection);
    gt = int8(gt);
end

if ~strcmpi(class(detection), class(gt))
    detection = int8(detection);
    gt = int8(gt);
end

if isa(detection, 'uint8') || isa(detection, 'uint16') || isa(detection, 'uint32') || isa(detection, 'uint64')
    detection = int8(detection);
end

if isa(gt, 'uint8') || isa(gt, 'uint16') || isa(gt, 'uint32') || isa(gt, 'uint64')
    gt = int8(gt);
end

if use_integral
    v = ver;
    if ~any(strcmp('Symbolic Math Toolbox', {v.Name}))
        use_integral = 0;
        ST = dbstack;
        warning(['Symbolic Math Toolbox not found, turning off numerical integral evaluation (results may not be as accurate as with this toolbox). To remove this warning turn off use_integral in ' ST(1).file]);
    end
    clear v ST
end



pos_label           = 1;
neg_label           = 0;
detection(detection == class_labels(2)) = pos_label;
gt(gt == class_labels(2))               = pos_label;
detection(detection == class_labels(1)) = neg_label;
gt(gt == class_labels(1))               = neg_label;




pi = sum(gt == pos_label) / numel(gt);


difference = gt - (detection*2);



tp = sum(difference == -positive_value);
fp = sum(difference == -(2*positive_value));
tn = sum(difference == neg_label);
fn = sum(difference == positive_value);




if (tp+fn) ~= 0
    tpr = tp / (tp+fn);
else
    tpr = 1;
end

if (fp+tn) ~= 0
    fpr = fp / (fp+tn);
else
    fpr = 0;
end

if (tp+fp) ~= 0
    p = tp / (tp+fp);
else
    p = 0;
end


r = tpr;

if tp ~= 0
    
    if use_integral
        
        % Integral -- slower but (slightly) more accurate
        
        q = @(t)(pi_function(t) ./ (pi_function(t) + ((1 - pi_function(t)) .* (pi ./ (1 - pi)) .* ((1 ./ p) - 1)))); % Eq. (4)
        p = (1/T) * integral(q, 0, T); % Eq. (4)
        
    else
        
        % Discrete -- faster but (slightly) less accurate
        
        t = 0 : T / 10 : T;
        p = (pi_function(t) ./ (pi_function(t) + ((1 - pi_function(t)) .* (pi ./ (1-pi)) .* ((1 ./ p) - 1)))); % Eq. (4)
        p = sum(p(~isnan(p))) ./ numel(p(~isnan(p)));       % Eq. (4)
        
    end
    
else
    
    p = 0;
    
end

end