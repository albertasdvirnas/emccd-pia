
function [regSizeThreshBlack , regSizeThreshWhite ] =  ...
        calc_region_size_threshold_using_clustering( regSizes , isRegBlack , allowedGapLength)

    %
    % Calculates region size thresholds for black and white
    % regions.
    %
    % Input:
    %
    % regionSizes = number pixels in each of the regions 
    % isRegBlack = parameter which is = 1 if a region is black, 
    %              and = 0 is white
    % allowedGapLength = allowed gap length. In sorted list of 
    %                    region sizes this, this gives us how large gaps, 
    %                    starting at region size = 1, we allow for the 
    %                    regions to be said to belong to 
    %                    the "background" cluster. Recommended value
    %                    is be of the order 1.
    %
    % Output:
    %
    % regSizeThreshBlack = region size threshold for black regions
    % regSizeThreshWhite = region size threshold for white regions
    %
    % Dependencies: None.
    %
  
                             
    %
    % Region size threshold for black regions 
   
    regSizesBlack = regSizes(find(isRegBlack));
    regSizesBlackUnique = unique(regSizesBlack); % returns a sorted array
    
    if regSizesBlackUnique(1) ~= 1   % if there are no regions of size = 1 at all, 
                                     % then we say that there are no background 
                                     % regions present whatsoever
          regSizeThreshBlack = 0;
     else
           
         % Go through the list of region sizes
         if length(regSizesBlackUnique) >= 2
             k = 1;
             gapLength = 0;
             while gapLength <= allowedGapLength & k < length(regSizesBlackUnique)
                 gapLength = regSizesBlackUnique(k+1) - regSizesBlackUnique(k); 
                 k = k+1;
             end
             regSizeThreshBlack =  regSizesBlackUnique(k-1);
         else
             regSizeThreshBlack = 1;
         end
         
    end
       
    
    %
    % Region size threshold for white regions 
    %  
    regSizesWhite = regSizes(find(~isRegBlack));    
    regSizesWhiteUnique = unique(regSizesWhite); % returns a sorted array
   
    if regSizesWhiteUnique(1) ~= 1   % if there are no regions of size = 1 at all, 
                                     % then we say that there are no background 
                                     % regions present whatsoever
          regSizeThreshWhite = 0;
     else
           
         % Go through the list of region sizes
         if length(regSizesWhiteUnique) >= 2
             k = 1;
             gapLength = 0;
             while gapLength <= allowedGapLength & k < length(regSizesWhiteUnique)
                 gapLength = regSizesWhiteUnique(k+1) - regSizesWhiteUnique(k); 
                 k = k+1;
             end
             regSizeThreshWhite =  regSizesWhiteUnique(k-1);
         else
             regSizeThreshWhite = 1;
         end
         
    end
     
  
         
end
