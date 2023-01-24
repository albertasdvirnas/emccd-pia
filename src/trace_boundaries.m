
function [boundariesCellArray , isBoundary] = trace_boundaries( labelIm  )

    %
    % Generates boundary curves (both exterior and interior) for all
    % regions from an region label matrix.
    %
    % Input:
    %
    % labelIm = matrix of the same size as the image where values
    %           give region labels
    % isBoundary = matrix of the same size as labelIm which contain a '1'
    %              if this pixel is a boundary pixel, else a '0'. Points at
    %              the image edges are here excluded.
    %
    % Output:
    %
    % boundariesCellArray = 
    %           cell array of cell arrays boundariesReg. 
    %           boundariesReg{k} gives a list of coordinates for boundary k
    %           of a given region (Note that a given region can have several
    %           disconnected boundaries).
    %       
    % Dependencies: None.
    %
    %


    % Hard-coded variable
    conn = 4; % connectivity


    % Pre-processing   
    nReg = max(labelIm(:));    % number of regions  
    [nRows , nCols ] = size(labelIm);
   
    % Find all boundary points all store in unordered lists   
    boundariesUnordered = cell(nReg,1);
    isBoundary = zeros(size(labelIm));
    for regIdx = 1:nReg  % loop over regions
        
         % Make an binarized image only containing 1s 
         % at positions corresponding to current region 
        binarizedImageReg = zeros(size(labelIm));
        idx = find(labelIm == regIdx);     
        binarizedImageReg(idx) = 1;
      
        boundaryList = zeros(numel(idx),2);    
        counter = 0;
        for k = 1:numel(idx) % loop over all pixels within a region
    
           [row,col] = ind2sub(size(labelIm),idx(k));
               
           isBoundaryImage = false;
           isBoundaryInterior = false;
           if (row == 1 | row == nRows | col == 1 | col == nCols)
                  % an image edge counts as a boundary
               isBoundaryImage = 1; 
           else     % find out if at least one neighbour is a '0' 
               nearestNeighbourBox = binarizedImageReg(row-1:row+1,col-1:col+1);
               if min(nearestNeighbourBox(:)) == 0
                   isBoundaryInterior = 1; 
               end
           end
           
           if isBoundaryImage | isBoundaryInterior 
               counter = counter + 1;
               boundaryList(counter,1) = row;
               boundaryList(counter,2) = col; 
              
           end
           
           if isBoundaryInterior 
               isBoundary(row,col) = 1;
           end
          
           
        end
        
        % Store result in a cell array
        boundariesUnordered{regIdx} = boundaryList(1:counter,:);   
                 
    end
 
    
    % Bring order to the boundary points   
    boundariesCellArray = cell(nReg,1);
   
    for regIdx = 1:nReg  % loop over regions
        
        % Make an binarized image only containing region regIdx pixels
        binarizedImageReg = zeros(size(labelIm));
        idx = find(labelIm == regIdx);     
        binarizedImageReg(idx) = 1;
        
        % Cell array where we store all boundaries for a given region
        boundariesReg = cell(1,1); 
         
         % A region can have many disconnected boundaries, let us find all
         % of these boundary curves
         boundaryPoints = boundariesUnordered{regIdx}; 
         boundaryPointsLinInd = sub2ind(size(labelIm),...
                 boundaryPoints(:,1),boundaryPoints(:,2) );
                           % convert boundary points to linear indices
         noOfBoundaryPoints = length(boundaryPointsLinInd);
          
         boundaryPointsExhaustedLinInd = []; 
         boundaryCounter = 1;
        
         while (length(boundaryPointsExhaustedLinInd) <  noOfBoundaryPoints) 
             
             [rows,cols] = ind2sub(size(labelIm),boundaryPointsLinInd);
           
             % Let us find a start point for the boundary tracing
             startRow = rows(1);
             startCol = cols(1);
             
             % Now trace out the boundary starting at (startRow,startCol).
             contour = bwtraceboundary(binarizedImageReg,[startRow , startCol],'E',conn);
             contourLinInd =  sub2ind(size(labelIm),contour(:,1),contour(:,2) );
             
             % Save contour
             boundariesReg{boundaryCounter} = contour;
                  
             % Remove the boundary points just detected from the list
             % to make sure we do not encounter this contour again
             boundaryPointsLinInd = setdiff(boundaryPointsLinInd,contourLinInd);
             
             % Keep track of all exhausted boundary points     
             boundaryPointsExhaustedLinInd = union(boundaryPointsExhaustedLinInd,contourLinInd);
           
             % Update counter
             boundaryCounter = boundaryCounter + 1;
          
         end
         boundariesCellArray{regIdx} = boundariesReg;
        
        

    end
    
   
    
    
end
