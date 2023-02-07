function [ binarizedImage ] = binarize_image_otsu( imageInput  )

    %
    % Binarizes an image using the Otsu method
    %
    % Input:
    %
    % imageInput = image
    %   
    % Output:
    %
    % binarizedImage = binarized image
    %
    % Dependencies: otsu/optimal_dim_bright_thresh_otsu.m
    % 
    
    % Find intensity threshold   
    [ intThresh , ~ ] = ...
            optimal_black_white_thresh_otsu(imageInput);
    disp(['Intensity thresh, Otsu = ',num2str(intThresh)])
        
    % Binarize
    binarizedImage = imbinarize(imageInput,intThresh);
  
   
end

