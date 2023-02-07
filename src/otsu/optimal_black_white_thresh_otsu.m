
function  [ intThresh ,idxThresh ] = ...
              optimal_black_white_thresh_otsu(imageInput)

    % Estimation an optimal intensity (image count) threshold,
    % using the Otsu method (minimize the sum of within-class variances)
    %
    % Input:
    %
    % imageInput = image (matrix) with intensity values
    %
    % Output:
    %
    % intThresh = intensity threshold value
    % idxThresh = index which corresponds to the threshold intensity
    %             in a list of sorted intensities,
    %             
    % 
    % Dependencies: None
    %
    
    % Otsu method for determining intensity threshold 
    imageMax = max(imageInput(:));
    normImage = imageInput/imageMax;
    normThreshold = graythresh(normImage);
    intThresh= normThreshold*imageMax;
  
    % Convert intensity threshold to index 
    % in a sorted list of  intensities 
    intVec = imageInput(:);
    intVecSorted = sort(intVec);  
    idx = find (intVecSorted <= intThresh);
    idxThresh = max(idx); 
    if idxThresh < 1   % some cleanup 
        idxThresh = 1;
    end
    if idxThresh > length(intVecSorted)
        idxThresh = length(intVecSorted);
    end
        
   
end
 