function detector_analysis(gts, annotations, detections, training_masks, gt_labels, detector_labels, evaluation_metric, pi_1, pi_2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETECTOR_ANALYSIS(GTS, DETECTIONS, TRAINING_MASKS, GT_LABELS, ...
%                   DETECTOR_LABELS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a set of ground truths, a set of detections, and
% (optionally) training masks depicting which part of the image was used
% for training the detectors.
%
% The detections are evaluated according to each of the ground truths, and
% then the detectors' ranks are calculated.
%
% If P-R curves are used then to remove the effects of the skew changes 
% between different ground truths the P-R curves are transformed to the 
% precision calculated at the mean skew of all the ground truths. For the 
% other methods the effects of skew are inherently normalised.
%
% See our paper for details regarding this methodology.
% 
% Input:
%
%      gts - a 3 dimensional matrix, n x m x f, where n and m are the
%            image's size and f is the number of ground truths. Within
%            which a one represents an object location and zero none.
%
%      detections - a 3 dimensional matrix, n x m x k, where n and m are
%                   the image's size and k is the number of detections.
%                   Within which a one represents a detection and zero
%                   none.
%
%      annotations - a 3 dimensional matrix, n x m x p, where n and m are
%                    the image's size and p is the number of annotations. 
%                    Within which a one represents an object's location and 
%                    zero none. These will be used as GTs, unless the
%                    matrix is empty.
%
%      training_masks - optional, a 3 dimensional integer matrix, 
%                       n x m x k, where n and m are the image's size and k
%                       is the number of detections. Within which a one 
%                       represents that the pixel was used for training and
%                       zeros that it was not (these positions are excluded
%                       from the performance evaluation). Ignored when
%                       using the Berkeley evaluation framework.
%
%      gt_labels - optional, a 1 x f cell vector that contains string 
%                  labels for each method used to calculate the ground 
%                  truths. If not GT labels are provided, default labels of
%                  GT1, ..., GTf will be generated.
%
%      detector_labels - optional, a 1 x k cell vector that contains string 
%                        labels for each detector used to calculate the 
%                        detections. If not detector labels are provided, 
%                        default labels of D1, ..., Dk will be generated.
%
%      evaluation_metric - specifies which metric to use (see SET VARIABLES
%                          section below) default is 'pr'.
%
%
% Output:
%
%      The function outputs the numerical results to the command prompt and 
%      also writes a file containing these results. The default name of the
%      output file is 'detector_ranks.txt' (which can be changed below).
%      The curve calculated according to each ground truth are saved as eps
%      and Matlab figure files. These will be located within a subdirectory 
%      called 'output' (which can be changed below) under the working 
%      directory.
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_filename    = 'detector_ranks.txt'; % Output filename for statistics

output_directory   = './output'; % Directory into which all the figures and 
                                 % text will be saved

if ~exist('evaluation_metric', 'var') % if no variable is pass then default
                                      % is used
    evaluation_metric  = 'pr';  % Which evaluation metric to use, options 
                                % are 'ipr', 'pr', 'roc', and 'berkeley'. 
                                % 'ipr' is best, but choose a reasonable 
                                % integrationrange, pi_1 and pi_2, below. 
                                % See the paper for an explanation.

    pi_1               = 0.01;
    pi_2               = 0.10;
    
end

interpolate        = 1;      % Interpole between points on the performance 
                             % curves, 0 for no or 1 for yes. Integration
                             % is slower but provides more accurate
                             % results.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path(path, genpath([fileparts(which('detector_analysis')) filesep 'functions']));

output_directory = correct_path(output_directory);
if ~isempty(output_directory) && ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

if ~exist('detector_labels', 'var')
    detector_labels = cell(1, size(detections, 3));
    for i = 1:numel(detector_labels); detector_labels{i} = ['D' num2str(i)]; end
end

if ~exist('training_masks', 'var')
	training_masks = false(size(detections));
end

if ~exist('gt_labels', 'var')
    gt_labels = cell(1, size(gts, 3));
    for i = 1:size(gts, 3); gt_labels{i} = ['GT' num2str(i)]; end
end

fid = fopen([output_directory output_filename], 'w');
finishup = onCleanup(@() fclose(fid));

gts(gts > 0) = 1;
annotations(annotations > 0) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THE PERFORMANCE OF THE DETECTORS ACCORING TO THE DIFFERENT 
% GROUND TRUTHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colours = lines(size(gts, 3)+2);

if ~isempty(annotations)
    % find the outliers (if any)
    [~, outliers] = cluster_annotators(annotations);
    
    % add the annotation labels to the list of GT labels
    for i = 1:size(annotations, 3)
        gt_labels = cat(2, gt_labels, ['A' num2str(i)]);
    end
end


if strcmpi(evaluation_metric, 'pr')
    skew = zeros(1, size(gts, 3)+size(annotations, 3));
    for i = 1:size(gts, 3)
        skew(i) = sum(sum(gts(:,:,i))) / numel(gts(:,:,i));
    end
    for i = 1:size(annotations, 3)
        skew(size(gts, 3)+i) = sum(sum(annotations(:,:,i))) / numel(annotations(:,:,i));
    end
    skew_m = mean(skew(:));
end


auc = zeros(size(detections, 3), size(gts, 3)+size(annotations, 3));
for j = 1:numel(detector_labels)
    
    h1 = figure('DefaultTextInterpreter', 'LaTex');
    a1 = axes('Parent', h1);
    
    curr_detection = detections(:,:,j);
    
    % Evaluation the detector's performance against each of the GTs
    for k = 1:size(gts, 3)

        curr_gt = gts(:,:,k);

        switch evaluation_metric
            case 'pr'
                skew = sum(curr_gt(:)) / numel(curr_gt(:));
                [ys, xs, auc(j,k)] = pr_curve(curr_detection(:), curr_gt(:), interpolate);
                [ys] = p_skew_transform(ys, skew, skew_m);
            case 'ipr'
                [ys, xs, auc(j,k)] = ipr_curve(curr_detection(~training_masks(:,:,j)), curr_gt(~training_masks(:,:,j)), interpolate, pi_1, pi_2);
            case 'roc'
                [ys, xs, auc(j,k)] = roc_curve(curr_detection(~training_masks(:,:,j)), curr_gt(~training_masks(:,:,j)), interpolate);
            case 'berkeley'
                [ys, xs, auc(j,k)] = ipr_curve_berkeley(curr_detection, curr_gt, interpolate, 0.0075, pi_1, pi_2, 1);
            case 'berkeley_internal'
                [ys, xs, auc(j,k)] = ipr_curve_berkeley_internal_only(curr_detection, curr_gt, interpolate, 0.0075, pi_1, pi_2, 1);
            otherwise
                error('Unknown evaluation metric');
        end
        
        plot(a1, xs, ys, 'Color', colours(k,:), 'LineWidth', 1.5);
        drawnow

        if k == 1
            hold on
        end
    
    end
    
    % Evaluate the detector's performance against the individual
    % annotations (if they exist)
    if ~isempty(annotations)
        for k = 1:size(annotations, 3)
            
            curr_gt = annotations(:,:,k);
            
            switch evaluation_metric
                case 'pr'
                    skew = sum(curr_gt(:)) / numel(curr_gt);
                    [ys, xs, auc(j,size(gts, 3)+k)] = pr_curve(curr_detection(:), curr_gt(:), interpolate);
                    [ys] = p_skew_transform(ys, skew, skew_m);
                case 'ipr'
                    [ys, xs, auc(j,size(gts, 3)+k)] = ipr_curve(curr_detection(~training_masks(:,:,j)), curr_gt(~training_masks(:,:,j)), interpolate, pi_1, pi_2);
                case 'roc'
                    [ys, xs, auc(j,size(gts, 3)+k)] = roc_curve(curr_detection(~training_masks(:,:,j)), curr_gt(~training_masks(:,:,j)), interpolate);
                case 'berkeley'
                    [ys, xs, auc(j,size(gts, 3)+k)] = ipr_curve_berkeley(curr_detection, curr_gt, interpolate, 0.0075, pi_1, pi_2, 1);
                case 'berkeley_internal'
                    [ys, xs, auc(j,size(gts, 3)+k)] = ipr_curve_berkeley(curr_detection, curr_gt, interpolate, 0.0075, pi_1, pi_2, 1);
                otherwise
                    error('Unknown evaluation metric');
            end
            
            if any(outliers == k)
                plot(a1, xs, ys, '--r', 'LineWidth', 1.5);
            else
                plot(a1, xs, ys, '--g', 'LineWidth', 1.5);
            end
            
            drawnow
            
            if k == 1
                hold on
            end
            
        end
    end
    

    axis(a1, 'square');
    axis(a1, [0 1 0 1]);
    title(a1, ['Evaluation of ' detector_labels{j} ' according to different GTs'], 'FontSize', 14);
    hold off
    legend(a1, strrep(gt_labels, '%', '\%'), 'Location', 'SouthEastOutside', 'Interpreter', 'latex');
    set(a1, 'FontSize', 14)
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
    drawnow
    saveas(h1, [output_directory 'detector_' detector_labels{j} '_gt_curves'], 'fig');
    saveas(h1, [output_directory 'detector_' detector_labels{j} '_gt_curves'], 'epsc');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THE RANKING OF THE DETECTORS MEASURE BY AUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(auc, 2)
    
    curr = auc(:,i);
    
    [~, order(:,i)] = sort(curr, 1, 'descend');
    
    fprintf('GT: %s\n', gt_labels{i});
    fprintf(fid, 'GT: %s\n', gt_labels{i});
    for j = 1:size(order, 1)
        fprintf('%s: %s (%f), ', num2str(j), detector_labels{order(j,i)}, curr(order(j,i)));
        fprintf(fid, '%s: %s (%f), ', num2str(j), detector_labels{order(j,i)}, curr(order(j,i)));
    end
    fprintf('\b\b\n\n');
    fprintf(fid, '\b\b\n\n');
end
fprintf('\b');
fprintf(fid, '\b');

[unique_orders, ~, ind] = unique(order', 'rows');
detect_order  = detector_labels(unique_orders);
for i = 1:size(detect_order, 1); col_headers{i} = ['Rank ' num2str(i)]; end

fprintf(fid, '\n');
fprintf('\n');
for i = sort(unique(ind))'

    gt_rank_names = gt_labels(ind == i);

    fprintf(fid, 'Rank %d: ', i);
    fprintf('Rank %d: ', i);
    for j = 1:numel(gt_rank_names)
        if j ~= numel(gt_rank_names)
            fprintf(fid, '%s, ', gt_rank_names{j});
            fprintf('%s, ', gt_rank_names{j});
        else
            fprintf(fid, '%s\n', gt_rank_names{j});
            fprintf('%s\n', gt_rank_names{j});
        end
    end

end

if strcmpi(evaluation_metric, 'ipr') || strcmpi(evaluation_metric, 'berkeley') || strcmpi(evaluation_metric, 'berkeley_internal')
    fprintf('\n(iPR integration range: pi_1 = %f, pi_2 = %f)\n', pi_1, pi_2);
    fprintf(fid, '\n(iPR integration range: pi_1 = %f, pi_2 = %f)', pi_1, pi_2);
end

delete(finishup);