function agreement_analysis(img, annotations, detections, training_masks, detector_labels, evaluation_metric, pi_1, pi_2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AGREEMENT_ANALYSIS(IMG, ANNOTATIONS, DETECTIONS, TRAINING_MASKS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes an image, a set of annotations marking objects within
% the image, (optionally) a set of detections calculated upon the image, 
% and (optionally) training masks depicting which part of the image was 
% used for training the detectors.
%
% Several statistics are calculated using this information:
% 
%   - Smyth's minimum bound on error;
%   - basic image feature and agreement correlations;
%   - agreement level percentage as a function of the number of annotators;
%   - outlier annotators;
%   - the sensitivity, specificity, positive predictive value, negative
%     predictive value, and kappa relative to 50% agreement between the
%     annotators;
%   - if detections are included, detector performance according to
%     agreement levels;
%   - if detections are included, correlations between detector output and
%     annotator agreement.
%
% See our paper for details regarding each of these methods.
% 
% Input:
%
%      img - the, n x m x b, image within which annotation is performed.
%            The image may have any number of bands, b. If the image has
%            three or more bands (b >= 3) then it is assumed that the first
%            three are the RGB bands. This can be changed by editing the 
%            rgb_bands variable (see below).
%
%      annotations - a 3 dimensional matrix, n x m x p, where n and m are
%                    the image's size and p is the number of annotations. 
%                    Within which a one represents an object's location and 
%                    zero none.
%
%      detections - optional, a 3 dimensional matrix, n x m x k, where n
%                   and m are the image's size and k is the number of 
%                   detections.Within which a one represents an object's 
%                   location and zero none.
%
%      training_masks - optional, a 3 dimensional integer matrix, 
%                   n x m x k, where n and m are the image's size and k is 
%                   the number of detections. Within which a one represents 
%                   that the pixel was used for training and zeros that it 
%                   was not (these positions are excluded from the 
%                   performance evaluation).
%
%      detector_labels - optional, a 1 x k cell vector that contains string 
%                  labels for each detector used to calculate the 
%                  detections. If not detector labels are provided, default 
%                  labels of D1, ..., Dk will be generated.
%
%      evaluation_metric - specifies which metric to use (see SET VARIABLES
%                          section below) default is 'pr'.
%
%
% Output:
%
%      The function outputs the numerical results to the command prompt and 
%      also writes a file containing these results. The default name of the
%      output file is 'agreement_statistics.txt' (which can be changed
%      below). The agreement curve plots are saved as eps and Matlab figure
%      files. These will be located within a subdirectory called 'output' 
%      (which can be changed below) under the working directory.
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_filename    = 'agreement_statistics.txt'; % Output filename for 
                                                 % statistics

output_directory   = './output'; % Directory into which all the figures and 
                                 % text will be saved


rgb_bands          = [1 2 3]; % If the image has three or more bands then
                              % this variable indicates which are the RGB 
                              % bands (in that order). This only concerns
                              % the calculation of the contrast image for
                              % correlation with agreement.

if ~exist('evaluation_metric', 'var') % if no variable is pass then default
                                      % is used

    evaluation_metric  = 'pr';  % Which evaluation metric to use, options 
                                % are 'ipr', 'pr', 'roc', and 'berkeley'. 
                                % 'ipr' is best, but choose a reasonable 
                                % integration range, pi_1 and pi_2, below. 
                                % See the paper for an explanation.

    pi_1               = 0.01;
    pi_2               = 0.10;

end

interpolate        = 0;      % Interpolate between points on the 
                             % performance curves, 0 for no or 1 for yes. 
                             % Integration is slower but provides more 
                             % accurate results.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




path(path, genpath([fileparts(which('agreement_analysis')) filesep 'functions']));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRE-PROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annotations(annotations > 0) = 1;
annotations = double(annotations);
agreement = sum(annotations, 3);

output_directory = correct_path(output_directory);
if ~isempty(output_directory) && ~exist(output_directory, 'dir')
    mkdir(output_directory);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE AGREEMENT MINIMUM BOUNDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_bound = smyth_min_error_bound(annotations);
fprintf('\nOverall errorbound according to Smyth: %f\n', error_bound);
fid = fopen([output_directory output_filename], 'w');
finishup = onCleanup(@() fclose(fid));
fprintf(fid, '\nOverall errorbound according to Smyth: %f\n', error_bound);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE IMAGE FEATURE AND AGREEMENT CORRELATION COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the correlation using each band of the image...

fprintf('\n');
fprintf(fid, '\n');
for i = 1:size(img, 3)
    
    [corr_c, p] = image_feature_agreement_correlation(img(:,:,i), agreement);
    
    fprintf('Image band %d and agreement correlation: %f (p = %f)\n', i, corr_c, p);
    fprintf(fid, 'Image band %d and agreement correlation: %f (p = %f)\n', i, corr_c, p);

end


% If the image is RGB then calculate the grey scale version...

if size(img, 3) == 3
    
    img_bw = rgb2gray(img);
    
    [corr_c, p] = image_feature_agreement_correlation(img_bw, agreement);
    
    fprintf('Grayscale image and agreement correlation: %f (p = %f)\n', corr_c, p);
    fprintf(fid, 'Grayscale image and agreement correlation: %f (p = %f)\n', corr_c, p);
    
    clear img_bw
end




% Calculate contrast correlation coefficient. The GT is adjusted to take
% the maximum agreement within the same neighbourhood as the contrast is 
% calculated.


if size(img, 3) >= 3
    img_rgb = img(:, :, rgb_bands);
else
    img_rgb = img;
end

img_cont = local_contrast(img_rgb);

gt_cont = agreement;
for i = 2:size(agreement,1)-1
    for j = 2:size(agreement, 2)-1
        gt_cont(i,j) = max(max(agreement(i-1:i+1, j-1:j+1)));
    end
end

[corr_c, p] = image_feature_agreement_correlation(img_cont, gt_cont);

fprintf('Image contrast and agreement correlation: %f (p = %f)\n', corr_c, p);
fprintf(fid, 'Image contrast and agreement correlation: %f (p = %f)\n', corr_c, p);

clear img_cont gt_cont img_rgb



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE AGREEMENT LEVEL PERCENTAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[percentages] = agreement_percentage(annotations);

fprintf('\n# annotators\t%% agreement\n');
fprintf(fid, '\n# annotators\t%% agreement\n');
for i = 1:numel(percentages)
    fprintf('%d\t\t%f\n', i, percentages(i));
    fprintf(fid, '%d\t\t%f\n', i, percentages(i));
end

h1 = figure('DefaultTextInterpreter', 'LaTex');
a1 = subplot(2,1,1);
plot(a1, 1:numel(percentages), percentages, '-ob', 'LineWidth', 2);
title(a1, 'Agreement level as a function of the number of annotators', 'FontSize', 14);
axis(a1, [1 numel(percentages) 0 100]);
xlabel(a1, 'Number of Annotators', 'FontSize', 14);
ylabel(a1, 'Agreement', 'FontSize', 14);
set(a1, 'XTick', 1:numel(percentages));
set(a1, 'FontSize', 14);
drawnow
saveas(h1, [output_directory 'agreement_level'], 'fig');
saveas(h1, [output_directory 'agreement_level'], 'epsc');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLUSTER ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[clustering, outliers] = cluster_annotators(annotations);

if ~isempty(outliers)
    fprintf('\nOutliers are: '); fprintf('A%d ', outliers); fprintf('\n');
    fprintf(fid, '\nOutliers are: '); fprintf(fid, 'A%d ', outliers); fprintf(fid, '\n');
else
    fprintf('\nNo Outliers Found.\n');
    fprintf(fid, '\nNo Outliers Found.\n');
end


h1 = figure('DefaultTextInterpreter', 'LaTex');
for i = 1:size(annotations, 3); lbls{i} = ['A' num2str(i)]; end
a1 = dendrogram(clustering, 'Labels', lbls, 'ColorThreshold', 0.9*max(clustering(:,3)));
for i = 1:numel(a1)
	set(a1(i), 'LineWidth', 2);
end
ylabel(gca, 'F$_1$-Score Difference', 'FontSize', 14);
xlabel(gca, 'Annotator', 'FontSize', 14);
title(gca, 'Dendrogram of annotator clustering', 'FontSize', 14);
set(gca, 'FontSize', 14, 'DefaultTextInterpreter', 'LaTex');
clear clustering
drawnow
saveas(h1, [output_directory 'dendrogram'], 'fig');
saveas(h1, [output_directory 'dendrogram'], 'epsc');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE SENSITIVITY, SPECIFICITY, PPV, NPV, AND KAPPA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[sens, spec, PPV, NPV, kappa] = annotator_statistics(annotations);


fprintf('\nAnn.\tSens.\t\tSpec.\t\tPPV\t\tNPV\t\tkappa\n');
fprintf(fid, '\nExp.\tSens.\t\tSpec.\t\tPPV\t\tNPV\t\tkappa\n');
for i = 1:size(annotations, 3)
    
    fprintf('%d\t%f\t%f\t%f\t%f\t%f\n', i, sens(i), spec(i), PPV(i), NPV(i), kappa(i));
    fprintf(fid, '%d\t%f\t%f\t%f\t%f\t%f\n', i, sens(i), spec(i), PPV(i), NPV(i), kappa(i));
    
end



if exist('detections', 'var')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE THE PERFORMANCE OF THE DETECTORS ACCORING TO DIFFERENT 
    % LEVELS OF AGREEMENT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('training_masks', 'var')
        training_masks = false(size(detections));
    end
    
    if ~exist('detector_labels', 'var')
        detector_labels = cell(1, size(detections, 3));
        for i = 1:numel(detector_labels); detector_labels{i} = ['D' num2str(i)]; end
    end

    fprintf('\nDetector\tCCO (p value) \t\tCCI (p value)\n');
    fprintf(fid, '\nDetector\tCCO\t\tCCI\n');
    CCO = zeros(numel(detector_labels), 1);
    CCI = zeros(numel(detector_labels), 1);
    pO  = zeros(numel(detector_labels), 1);
    pI  = zeros(numel(detector_labels), 1);
    for j = 1:numel(detector_labels)
        
        curr_detection = detections(:,:,j);
        
        switch evaluation_metric
            case 'pr'
                [~,~,~,~,h1,a1] = agreement_curves(curr_detection(~training_masks(:,:,j)), agreement(~training_masks(:,:,j)), interpolate, evaluation_metric);
            case 'ipr'
                [~,~,~,~,h1,a1] = agreement_curves(curr_detection(~training_masks(:,:,j)), agreement(~training_masks(:,:,j)), interpolate, evaluation_metric, pi_1, pi_2);
            case 'roc'
                [~,~,~,~,h1,a1] = agreement_curves(curr_detection(~training_masks(:,:,j)), agreement(~training_masks(:,:,j)), interpolate, evaluation_metric);
            case 'berkeley'
                [~,~,~,~,h1,a1] = agreement_curves(curr_detection(~training_masks(:,:,j)), agreement(~training_masks(:,:,j)), interpolate, evaluation_metric, pi_1, pi_2);
            case 'berkeley_internal'
                [~,~,~,~,h1,a1] = agreement_curves(curr_detection(~training_masks(:,:,j)), agreement(~training_masks(:,:,j)), interpolate, evaluation_metric, pi_1, pi_2);
        end
        
        l1 = legend;
        set(l1, 'Location', 'SouthEastOutside');
        title(a1, ['Agreement Curves of ' detector_labels{j}], 'FontSize', 14);
        drawnow
        saveas(h1, [output_directory 'detector_' detector_labels{j} '_agreement_curves'], 'fig');
        saveas(h1, [output_directory 'detector_' detector_labels{j} '_agreement_curves'], 'epsc');
        
        [CCO(j), pO(j)] = corr(double(curr_detection(agreement > 0 & ~training_masks(:,:,j))), double(agreement(agreement > 0 & ~training_masks(:,:,j))));
        [CCI(j), pI(j)] = corr(double(curr_detection(~training_masks(:,:,j))), double(agreement(~training_masks(:,:,j))));

        fprintf('%s\t\t%f (%f)\t%f (%f)\n', detector_labels{j}, CCO(j), pO(j), CCI(j), pI(j));
        fprintf(fid, '%s\t\t%f (%f)\t%f (%f)\n', detector_labels{j}, CCO(j), pO(j), CCI(j), pI(j));

    end

    clear curr_detection
    
end

if strcmpi(evaluation_metric, 'ipr') || strcmpi(evaluation_metric, 'berkeley')
    fprintf('\n(iPR integration range: pi_1 = %f, pi_2 = %f)\n', pi_1, pi_2);
    fprintf(fid, '\n(iPR integration range: pi_1 = %f, pi_2 = %f)', pi_1, pi_2);
end


delete(finishup);

end