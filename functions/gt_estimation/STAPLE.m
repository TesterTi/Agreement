function staple_gt = STAPLE(annotations)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SENS, SPEC, PPV, NPV, KAPPA] = ANNOTATOR_STATISTICS(ANNOTATIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a set of annotations and the ground truth according 
% to the STAPLE method.
%
% The STAPLE method is taken from:
%
%   B. Landman, J. Bogovic and J. Prince, 'Simultaneous truth and 
%        performance level estimation with incomplete, over-complete, and 
%        ancillary data', in: Proc. SPIE Medical Imaging 2010: Image 
%        Processing, vol. 7623, 2010.
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
%      staple_gt - a n x m matrix representing the resulting ground truth. 
%                  Within which a one represents an object's location and 
%                  zero none.
%
%
%	T. Lampert, A. Stumpf, and P. Gancarski, 'An Empirical Study into 
%       Annotator Agreement, Ground Truth Estimation, and Algorithm 
%       Evaluation', IEEE Transactions on Image Processing 25 (6): 
%       2557–2572, 2016.
%
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    memory_allocation = '2048'; % Memory allocation to the java virtual 
                                % machine for ImageJ, if you receive some 
                                % out of memory errors increase this, if 
                                % you receive some error about not being 
                                % able to allocate memory, then reduce this.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                

    annotations(annotations > 0) = 1;
    annotations = uint8(annotations);
    
    imj_path     = correct_path([fileparts(which('STAPLE')) filesep 'support_software/ImageJ']);
    staple_path  = correct_path([fileparts(which('STAPLE')) filesep 'support_software/STAPLE']);
    working_path = correct_path([fileparts(which('STAPLE')) filesep 'support_software']);

    switch computer
        case 'MACI64'
            staple_path = [staple_path 'OSX' filesep];
        case 'PCWIN'
            staple_path = [staple_path 'Windows' filesep];
        case 'PCWIN64'
            staple_path = [staple_path 'Windows' filesep];
        case 'GLNXA64'
            staple_path = [staple_path 'OSX' filesep];
        case 'GLNX86'
            staple_path = [staple_path 'OSX' filesep];
    end

    if system('java -version') ~= 0
        error('Java not found! The Java virtual machine is not installed or the system''s path is not setup correctly.');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONVERT IMAGES TO NRRD FORMAT USING IMAGEJ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    imj_macro = [];
    for i = 1:size(annotations,3)

        imwrite(annotations(:,:,i), [working_path 'temp' num2str(i) '.tif'], 'TIF');

        imj_macro = [imj_macro 'path = "' strrep(working_path, filesep, [filesep filesep]) 'temp' num2str(i) '.tif"; open(path); run("Nrrd Writer", "nrrd=' strrep(working_path, filesep, [filesep filesep]) 'temp' num2str(i) '.nrrd"); close();'];

    end
    imj_macro = [imj_macro 'run(''Quit'');'];

    % Write macro to file so that we can execute ImageJ without a GUI
    fid = fopen([working_path 'temp.macro'], 'w');
    finishup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', imj_macro);
    delete(finishup);

    imj_cmd = ['java -jar -Xmx' memory_allocation 'm ' imj_path 'ij.jar -ijpathpat ' imj_path 'plugins -batchpath "' working_path 'temp.macro"'];

    if system(imj_cmd) ~= 0
        error('An error occurrect while executing ImageJ. Please check the above error message.');
    end


    for i = 1:size(annotations,3)
        % ImageJ always puts 3D co-ordinates for a 2D image so need to correct
        % that

        fid = fopen([working_path 'temp' num2str(i) '.nrrd']);
        finishup = onCleanup(@() fclose(fid));
        A = fread(fid, '*char');
        delete(finishup);

        B = strrep(A', '(1.0,0,0) (0,1.0,0)', '(1.0,0) (0,1.0)');
        fid = fopen([working_path 'temp' num2str(i) '_new.nrrd'], 'w');
        finishup = onCleanup(@() fclose(fid));
        fwrite(fid, char(B)');
        delete(finishup);

        movefile([working_path 'temp' num2str(i) '_new.nrrd'], [working_path 'temp' num2str(i) '.nrrd']);

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE STAPLE GT USING STAPLE EXECUTABLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    staple_cmd = [staple_path 'crlSTAPLE --assignConsensusVoxels 0 -o ' working_path 'weights.nrrd '];
    for i = 1:size(annotations,3)
        staple_cmd = [staple_cmd working_path 'temp' num2str(i) '.nrrd '];
    end

    if system(staple_cmd) ~= 0
        error('An error occurrect while executing STAPLE. Please check the above error message.');
    end

    % Need to convert first output into segementation

    staple_cmd = [staple_path 'crlIndexOfMaxComponent ' working_path 'weights.nrrd ' working_path 'GT_staple.nrrd'];

    if system(staple_cmd) ~= 0
        error('An error occurrect while executing STAPLE. Please check the above error message.');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONVERT BACK TO IMAGE FORMAT USING IMAGEJ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    imj_macro = ['run("Nrrd Reader", "load=' strrep(working_path, filesep, [filesep filesep]) 'GT_staple.nrrd"); saveAs("Tiff", "' strrep(working_path, filesep, [filesep filesep]) 'GT_staple.tif"); run("Quit");'];

    % Write macro to file so that we can execute ImageJ without a GUI
    fid = fopen([working_path 'temp.macro'], 'w');
    finishup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', imj_macro);
    delete(finishup);

    imj_cmd = ['java -jar -Xmx' memory_allocation 'm ' imj_path 'ij.jar -ijpathpat ' imj_path 'plugins -batchpath "' working_path 'temp.macro"'];

    if system(imj_cmd) ~= 0
        error('An error occurrect while executing ImageJ. Please check the above error message.');
    end



    %%%%%%%%%%%%%%%%%%%%%%%
    % Read in STAPLE result
    %%%%%%%%%%%%%%%%%%%%%%%

    staple_gt = imread([working_path 'GT_staple.tif']);


    delete([working_path 'GT_staple.nrrd']);
    delete([working_path 'GT_staple.tif']);
    delete([working_path 'weights.nrrd']);
    delete([working_path 'temp.macro']);
    for i = 1:size(annotations,3)
        delete([working_path 'temp' num2str(i) '.nrrd']);
        delete([working_path 'temp' num2str(i) '.tif']);
    end

end