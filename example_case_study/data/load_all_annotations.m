function annotations = load_all_annotations


for i = 1:13
    
    fprintf('%3d of 13', i);
    
    x = imread(['gt_A' num2str(i) '.tif']);
    annotations(:,:,i) = x;
    
    fprintf('\b\b\b\b\b\b\b\b\b');
end