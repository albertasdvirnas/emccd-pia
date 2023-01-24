function [mergedLabelIm , origLabels] = ...
    generate_final_segmented_regions(labelIm, isRegBlack, regSizeThreshBlack , regSizeThreshWhite )

    %
    % Generates final segmented regions.
    %
    % Input:
    %
    % labelIm = matrix of the same size as the image where values
    %           give region labels 
    % isRegBlack = parameter which if = 1 if a region is black, else = 0.
    % regSizeThreshBlack = region size threshold for black regions 
    %                       within white regions 
    % regSizeThreshWhite = region size threshold for white regions 
    %                      within black regions 
    %         
    % Output:
    %
    % mergedLabelIm = each pixel is assigned a region label stored in this
    %                 matrix   
    % origLabels =  cell array which tells which set of original regions 
    %               which are contained in each new region.
    %               Example: origLabels{3} = [4 8], tells 
    %               that the merged region 3 is composed of
    %               the regions 4 and 8 in the original (input) image. 
    %
    % Dependencies: None.
    %
    
    % Hard-coded variable:
    conn = 4; % connectivity
      
    % Number of regions (both white and black regions) 
    % in original image
    nRegOriginal = max(labelIm(:));
 
    
    % Binarized image where we retained all white regions
    % 
    
    blackRegIdx = find(isRegBlack == 1);
    binarizedImage = 1 - ismember(labelIm,blackRegIdx);
%    % Alternative: much slower
%     binarizedImage = zeros(size(labelIm));
%     for idxReg = 1:nRegOriginal
%     
%        if ~isRegBlack(idxReg)
%            idx = find(labelIm == idxReg);
%            binarizedImage(idx) = 1;
%        end
%        
%     end
%     
    
    
    
    % Find small connected components of white pixels within
    % black regions and remove these
    CC = bwconncomp(binarizedImage,conn);  
    nReg = CC.NumObjects;  % number of regions
    regSizes = cellfun(@numel,CC.PixelIdxList); % number of pixels in each regions
    for idxReg = 1:nReg
        if regSizes(idxReg) <= regSizeThreshWhite
            binarizedImage(CC.PixelIdxList{idxReg}) = 0;
        end
    end
  
   
    % Find  connected components of black pixels within
    % white regions and remove these 
    binarizedImageInv = 1 - binarizedImage;
    CC = bwconncomp(binarizedImageInv,conn);  
    nReg = CC.NumObjects;  % number of regions
    regSizes = cellfun(@numel,CC.PixelIdxList); % number of pixels in each regions
    for idxReg = 1:nReg
        if regSizes(idxReg) <= regSizeThreshBlack
            binarizedImageInv(CC.PixelIdxList{idxReg}) = 0;
        end
    end
    binarizedImage = 1 - binarizedImageInv;
  
      
    % Find connected components in final black-and-white image
    CC = bwconncomp(binarizedImage,conn);   
    mergedLabelIm = double(labelmatrix(CC));

    
    % Figure out which original regions each new region is composed of.
    nReg = max(mergedLabelIm(:)); % number of merged regions
    origLabels = cell(nReg,1);
    for k = 1:nReg
         idx = find(mergedLabelIm == k);
         origLabelsAll = labelIm(idx);
         origLabelsUnique = unique(origLabelsAll);
         origLabelsUnique = sort(origLabelsUnique); 
         origLabels{k} = origLabelsUnique;
    end
    
      
    
end
