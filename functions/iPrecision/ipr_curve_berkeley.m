function [p, r, au, threshold_values] = ipr_curve_berkeley(response, gt, interpolate, maxDist, pi1, pi2, thin_gt)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [p, r, au, threshold_values] = ipr_curve_berkeley(response, gt, ...
%               maxDist, interpolate, pi1, pi2)
%
% Calculate iPrecision/Recall curve using the Berkeley matching criteria.
% The variable thin_gt indicates whether the ground truth should be thinned
% which is useful when ground truth methods don't leave 1 pixel wide GT
% estimates (default is off). This requires the Image Processing Toolbox.
%
% INPUT
%	response    : An N x M Detector response
%	gt          : Ground Truth image (also N x M)
%   maxDist     : For corresponding pixels, (0.0075), see Berkeley
%                 evaluation documentation.
%   interpolate : interpolate between iPR points (0 or 1)
%   pi1         : iPR integration limit
%   pi2         : iPR integration limit
%
% OUTPUT
%	p                : iPrecision
%	r                : Recall
%	au               : Area under the curve
%	threshold_values : Threshold values used to derive the iPR curve
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
number_of_data_points = 100;       % of the PR Curve

if ~exist('thin_gt', 'var')
    thin_gt = 0;
end
 
use_integral = 1;       % calculate using integral or discrete (integral is 
                        % slower but more accurate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~isscalar(pi1) || ~isfloat(pi1) || pi1 > 1 || pi1 < 0
    error('Pi 1 must be a floating point scalar between 0 and 1');
end

if ~isscalar(pi2) || ~isfloat(pi2) || pi2 > 1 || pi2 < 0
    error('Pi 2 must be a floating point scalar between 0 and 1');
end

if isvector(response)
    error('The algorithm''s response should be in the format of a N x M matrix, where N is the image height and M is its width.');
end

if isvector(gt)
    error('The gt should be in the format of a N x M matrix, where N is the image height and M is its width.');
end

if all(size(response) ~= size(gt(:,:,1)))
    error('The detection and gt should be of the same size');
end

positive_value = 1;
gt(gt > 0)     = positive_value;

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

if thin_gt
    v = ver;
    if ~any(strcmp('Image Processing Toolbox', {v.Name}))
        thin_gt = 0;
        ST = dbstack;
        warning(['Image Processing Toolbox not found, turning off ground truth thinning (performance will not be as expected). To remove this warning turn off thing_gt in ' ST(1).file]);
    end
    clear v ST
end

berkeley_path     = correct_path([fileparts(which('ipr_curve_berkeley')) filesep 'support_functions']);
path(path, berkeley_path);




threshold_values = linspace(double(min(min(response))), double(max(max(response))), number_of_data_points);

% zero all counts
p  = zeros(1, number_of_data_points);
tp = zeros(1, number_of_data_points);
fp = zeros(1, number_of_data_points);
tn = zeros(1, number_of_data_points);
fn = zeros(1, number_of_data_points);
cntR = zeros(size(threshold_values));
sumR = zeros(size(threshold_values));

for i = 1:numel(threshold_values)
        
    detection = zeros(size(response), 'int8');

    detection(response >= threshold_values(i)) = 1;

    accP = zeros(size(response));
    
    for j = 1:size(gt, 3)
            
        curr_gt = gt(:, :, j);
            
        if thin_gt
            curr_gt = double(bwmorph(curr_gt, 'thin', inf));
        end

        [match1, match2] = correspondPixels(double(detection), double(curr_gt), maxDist);

        accP = accP | match1;

        % compute recall
        sumR(i) = sumR(i) + sum(curr_gt(:));
        cntR(i) = cntR(i) + sum(match2(:) > 0);

    end
        
	sumP = sum(detection(:));
        
    tn(i) = tn(i) + sum(sum(accP == 0 & curr_gt == 0));
    fn(i) = fn(i) + sum(sum(accP == 0 & curr_gt > 0));
    
    tp(i) = sum(match1(:));
    fp(i) = sumP - tp(i);
    
    phi = sum(gt(:) > 0) / sum(gt(:) == 0);                                % !!!NEED TO CHECK HOW TO CALCULATE DATASET SKEW!!!
    
    % compute the iPR point
    if use_integral
        
        % Integral -- slower but (slightly) more accurate
        q = @(pi)((pi * tp(i)) ./ ((pi * tp(i)) + ((1-pi) * phi * fp(i)))); 
        p(i) = (1/(pi2-pi1))*integral(q, pi1, pi2); % Eq. (7)
        
    else
        
        % Discrete -- faster but (slightly) less accurate
        pi = pi1 : (pi2 - pi1) / 1000 : pi2;
        p(i) = (pi * tp(i)) ./ ((pi * tp(i)) + ((1-pi) * phi * fp(i)));
        p(i) = sum(p(~isnan(p))) / numel(p(~isnan(p)));       % Eq. (7)
        
    end
    
	%[~, ~, p(i), ~, tp(i), fp(i), tn(i), fn(i)] = ipr_point(accP(:), curr_gt(:), pi1, pi2);
    
end

r = cntR ./ sumR;




% Interpolate along the curve if desired
if interpolate
    
    ps_n = [];
    rs_n = [];
    
    for j = 1:number_of_data_points-1
        
        [ps1, rs1] = ipr_interpolate(tp(j+1), tp(j), fp(j+1), fp(j), tn(j), fn(j), pi1, pi2);
        
        ps_n = [ps_n, p(j), ps1(end:-1:1)];
        rs_n = [rs_n, r(j), rs1(end:-1:1)];

    end
    
    p = ps_n;
    r = rs_n;
    
end



% Calculate the curve's AUC
au = auc(p(~isnan(p)), r(~isnan(p)));

end