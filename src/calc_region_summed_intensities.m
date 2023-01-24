
function summedInt = calc_region_summed_intensities( imageInput , labelIm  )

    %
    % Calculates the sum of intensities in each region
    %
    % Input:
    %
    % image =  image (or matrix in double precision) 
    % labelIm = matrix of the same size as the image where values
    %           give region labels
    %
    % Output:
    %
    % summedInt = sum of intensities over all pixels in a region.
    %
    % Dependencies: None.
    %
  
    nReg = max(labelIm(:));    % number of regions 
    summedInt = zeros(1,nReg); 
    for idx = 1:nReg 
        regIdx = (labelIm == idx);
        intensities = imageInput(regIdx);
        summedInt(idx) = sum(intensities);      
    end
       
end
