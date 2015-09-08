% Builds the support functions for the function ipr_curve_berkeley

%output_path     = [fileparts(which('build')) filesep 'toolbox' filesep 'support_functions' filesep];
output_path = ['.' filesep 'functions' filesep 'iPrecision' filesep 'support_functions' filesep];
source_path = [output_path 'source' filesep];
cd(source_path)
mex -largeArrayDims -v CXXFLAGS="\$CXXFLAGS -O3 -DNOBLAS" -outdir ../ correspondPixels.cc csa.cc kofn.cc match.cc Exception.cc Matrix.cc Random.cc String.cc Timer.cc
cd ../../../..

clear output_path source_path