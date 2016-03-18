function [xs, ys, aucs, thresholds, h1, a1] = agreement_curves(detection, agreement, interpolate, evaluation_metric, pi_1, pi_2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xs, ys, aucs, thresholds, h1, a1] = AGREEMENT_CURVES(DETECTION, GT, ... 
%         INTERPOLATE, EVALUATION_METRIC, PI_1, PI_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a detection and an agreement map and calculates the
% perfrmance curves of the detector (which resulted in the detection)
% occording to each level of agreement.
%
% If P-R curves are used then to remove the effects of the skew changes 
% between different ground truths (levels of agreement) the P-R curves are
% transformed to the precision calculated at the mean skew of all the
% agreement ground truths. For the other methods the effects of skew are
% inherently normalised.
%
% See our paper for details regarding each of these methods.
% 
% Input:
%
%      detection - a n x m matrix, where n and m are the image's size. 
%                  Within which a one represents the detection of an object
%                  and zero none.
%
%      agreement - a n x m matrix where n and m are the image's size. 
%                  Within which each pixel contains a count of the number 
%                  of times an annotator marked that position. 
%
%      interpolate - a binary scalar, 1 means to interpolate between points
%                    on the performance curve, or 0 not to.
%
%      evaluation_metric - 'ipr', 'pr', 'roc', or 'berkeley'. The type of 
%                          performance curve to be used:
%                           - ipr, integrated precision-recall
%                           - pr, precision-recall
%                           - roc, receiver operating characteristic
%                           - berkeley, ipr curve using the berekely 
%                             evaluation framework (see paper for details).
%
%      pi_1 and pi_2 - are only used if the type is 'ipr' or 'berkeley'.
%                      These specify the skew ranges within which to 
%                      integrate (see the paper for details)
%
%
% Output:
%
%      xs and ys - two N x M matrices, each row contains a vector of 
%                  coordinate points (x and y respecively) of the M curves 
%                  calculated. Where M is the maximum agreement found 
%                  within the agreement matrix.
%
%      aucs - a 1 x M vector containing the area under the curves
%             calculated.
%
%      thresholds - a 1 x N vector of the thresholds applied to the 
%                   detection to calculate each curve. Where N is the
%                   number of points on each curve.
%
%      h1 and a1 - handles to the figure and axis respectively.
%
% The function also outputs a figure representing the performance curves.
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




agreement = double(agreement);
detection = double(detection);



cc = jet(max(agreement(:)));

h1 = figure('DefaultTextInterpreter', 'LaTex');
a1 = axes('Parent', h1);

% If P-R curves are desired we need to remove the effects of the different
% skews encountered when creating the GTs and therefore we take the means
% skew and adjust all the P-R curves to that.
if strcmpi(evaluation_metric, 'pr')
    skew = zeros(1, max(agreement(:)));
    for i = 1:max(agreement(:))
        skew(i) = sum(agreement >= i) / numel(agreement);
    end
    skew_m = mean(skew(:));
end

ys = zeros(max(agreement(:)), 1);
xs = zeros(max(agreement(:)), 1);
thresholds = zeros(max(agreement(:)), 1);
aucs = zeros(max(agreement(:)), 1);
for i = 1:max(agreement(:))
    if any(agreement(:) == i)
        
        gt = zeros(size(agreement));
        gt(agreement >= i) = 1;
        
        switch evaluation_metric
            
            case 'pr'
                t_skew = sum(gt(:)) / numel(gt);
                [t_ys, t_xs, aucs(i), t_thresholds] = pr_curve(detection(:), gt(:), interpolate);
                [t_ys] = p_skew_transform(t_ys, t_skew, skew_m);
                
            case 'ipr'
                [t_ys, t_xs, aucs(i), t_thresholds] = ipr_curve(detection(:), gt(:), interpolate, pi_1, pi_2);
                
            case 'roc'
                [t_ys, t_xs, aucs(i), t_thresholds] = roc_curve(detection(:), gt(:), interpolate);
            
            case 'berkeley'
                [t_ys, t_xs, aucs(i), t_thresholds] = ipr_curve_berkeley(detection, gt, interpolate, 0.0075, pi_1, pi_2, 1);
                
            case 'berkeley_internal'
                [t_ys, t_xs, aucs(i), t_thresholds] = ipr_curve_berkeley_internal_only(detection, gt, interpolate, 0.0075, pi_1, pi_2, 1);
                
            otherwise
                error('Unknown evaluation metric!');
                
        end
        
        ys(i, 1:numel(t_ys)) = t_ys;
        xs(i, 1:numel(t_xs)) = t_xs;
        thresholds(1, 1:numel(t_thresholds)) = t_thresholds;
        
        plot(a1, xs(i,:), ys(i,:), 'color', cc(i,:), 'LineWidth', 2);
        
        if i == 1
            hold on
        end
        drawnow
        
    end
end


axis(a1, [0 1 0 1]);
axis(a1, 'square')
title(a1, evaluation_metric, 'FontSize', 14);
switch evaluation_metric
    
    case 'pr'
        xlabel(a1, 'Recall', 'FontSize', 14);
        ylabel(a1, 'Precision', 'FontSize', 14);
        
    case 'ipr'
        xlabel(a1, 'Recall', 'FontSize', 14);
        ylabel(a1, '$\bar{\textrm{P}}$recision', 'FontSize', 14);
    
    case 'roc'
        xlabel(a1, 'FPR', 'FontSize', 14);
        ylabel(a1, 'TPR', 'FontSize', 14);
    
    case 'berkeley'
        xlabel(a1, 'Recall', 'FontSize', 14);
        ylabel(a1, '$\bar{\textrm{P}}$recision', 'FontSize', 14);
    
end
set(a1, 'XTick', [0:0.2:1])
set(a1, 'YTick', [0:0.2:1])
set(a1, 'FontSize', 14)
hold off

legend_labels = cell(1, max(agreement(:)));
for i = 1:max(agreement(:))
    legend_labels{i} = ['$\ge ' num2str(i) '$'];
end
legend(a1, legend_labels, 'FontSize', 14, 'Interpreter', 'latex');

end