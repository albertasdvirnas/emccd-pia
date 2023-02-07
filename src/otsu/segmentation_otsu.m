function segOutputOtsu = segmentation_otsu(images)
   
    %
    % Perform image segmentation using the Otsu method
    % 
    % Input:
    % 
    % images = struct with images and associated information
    % 
    % Output:
    %
    % segOtsuOutput = struct containing all output of the segmentation
    %
    % Dependencies: otsu/binarize_image_otsu.m
    %
    
    % Hard-coded variables
    conn = 4;  % connectivity used when finding connected components
       
    % struct containing the segmentation output
    segOutputOtsu = struct();
    
    % Binarize
    [ binarizedImage ] = binarize_image_otsu( images.imAverage  );

    % Find white regions (connected components) in binarized image 
    CC = bwconncomp(binarizedImage,conn);   
    labelIm = double(labelmatrix(CC));
   
    % Return output
    segOutputOtsu.labelIm = labelIm;
    segOutputOtsu.binarizedImage = binarizedImage;
    
       
        
end
