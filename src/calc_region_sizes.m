
function [regSize] = calc_region_sizes( labelIm  )

    %
    % Calculates the number of pixels in each region
    %
    % Input:
    %
    % labelIm = matrix of the same size as the image where values
    %           give region labels
    %
    % Output:
    %
    % regSize = number pixels in each of the regions 
    %
    % Dependencies: None.
    %
  
    
    % Make a histogram
    nReg = max(labelIm(:));    % number of regions
    binEdges = 0.5:1:nReg+0.5;  % bin edges
    regSize = histcounts(labelIm(:),binEdges);
    

%     % Much slower variant: 
%     for idx = 1:nReg 
%         regIdx = find(labelIm == idx);
%         regSize(idx) = length(regIdx);            
%     end


end
