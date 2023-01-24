function [labelIm, isRegBlack ] = ...
             find_regions_in_binarized_image( binarizedImage )

    %
    % Uses a binarized image as input and groups connected
    % white pixels and connected black pixels respectively into
    % regions. Each region is assigned a unique label and is also 
    % assigned a 'black' or 'white' label.
    %
    % Input:
    %
    % binarizedImage = binarized image
    %
    % Output:
    %
    % labelIm = each pixel's label is stored in this matrix 
    %           (matrix size = size of image).
    %           Black regions are given labels 1,...,nRegBlack, and
    %           white regions are given labels nReg+1,...,max(labelIm(:)).
    %           where "region 1" is the largest connected component.
    % isRegBlack = a vector of the same length as the number of regions. 
    %             an element = 1 says that this region is black
    %             and an element = 0 says that this region is white.
    %
    % Dependencies: None.
    % 
     
    % Hard-coded variables
    conn = 4;  % connectivity used when finding connected components
       
    % White regions in binarized image 
    CC = bwconncomp(binarizedImage,conn);   
    labelImWhite = double(labelmatrix(CC));
    nRegWhite = CC.NumObjects;
    
    % Black regions in binarized image
    binarizedImageInv = 1 - binarizedImage;
    CC = bwconncomp(binarizedImageInv,conn);  %
    labelImBlack = double(labelmatrix(CC));
    nRegBlack = CC.NumObjects;
    
    % Swap labels so that the largest black connected component 
    % has region label = 1  
    idxPreviousRegionOne = find(labelImBlack == 1); 
                              % Store previous "region 1" pixels indices
    regSizesBlack = cellfun(@numel,CC.PixelIdxList); % all region sizes
    [~ ,idx] = max(regSizesBlack);  % find largest connected component
    labelImBlack(CC.PixelIdxList{idx}) = 1;
    labelImBlack(idxPreviousRegionOne) = idx;
    
    
    % Merge the white and black regions 
    labelIm = labelImBlack;
    idx = find(binarizedImage == 1);
    labelIm(idx) = labelImWhite(idx) + nRegBlack;
   
    % Store info about whether a given region is black or white
    isRegBlack = zeros(1,nRegBlack +  nRegWhite);
    isRegBlack(1:nRegBlack) = 1;
    
   
   
end

