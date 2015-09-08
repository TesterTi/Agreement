function p = correct_path(p)

% Makes sure that the path has the correct file separators for the OS and
% that it is terminated in a file separator.

p = strrep(p, '\', filesep);
p = strrep(p, '/', filesep);
if ~strcmp(p(end), filesep)
    p = [p, filesep];
end
