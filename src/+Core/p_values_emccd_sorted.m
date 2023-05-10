
function [pVals,pdfUnique,cdfsUnique,intUnique] = p_values_emccd_sorted(sortedInt, lambda, gain, adFactor, countOffset, roNoise)

    %
    % Calculates p-values using the EMCCD distribution
    % The input intensity values are assumed to be sorted.
    %
    % Input:
    %
    % intensities = vector with sorted intensity values
    % chipPars = struct containing the chip parameters
    %
    % Output:
    %
    % pVals = p-values for each input intensity
    %
    % Dependencies: emccd_distribution/pdf_cdf_emccd.m 
    % 
    
    % Hard-coded variable
%     N = 2^8;     % number of integration points when calculating 
                 % PDF through  numerical inverse Fourier-transform
                 % of the characteristic function   

    import Core.calc_bounds;
    [L,U,EX,STD] = calc_bounds(lambda,gain,adFactor,countOffset,roNoise);
  
   % Turn into a row vector
    sortedInt = sortedInt(:);  
    
    % Find all unique intensity values
   [intUnique,idx]  = unique(sortedInt);
   
   % take only those smaller than integration range
    idx = idx(intUnique <= floor(U));
    intUnique = intUnique(intUnique <= floor(U));
    
    nVals = length(idx);

    import Core.pdf_cdf_emccd;

    % Evaluate the CDF only at the unique intensities 
    if ~isempty(intUnique)
        [pdfUnique,cdfsUnique] = pdf_cdf_emccd(intUnique', lambda, gain, adFactor, countOffset, roNoise, L, U);
        % remove ones out of range
        cdfsUnique = min(cdfsUnique,1);
        cdfsUnique = max(cdfsUnique,0);

        % Calculate p-values for all input intensities
        pVals = zeros(1,length(sortedInt));
        for k=1:nVals-1   
            pVals(idx(k):idx(k+1)-1) = 1 - cdfsUnique(k);
        end
        pVals(idx(nVals):end) = 1 - cdfsUnique(end);
    else
        pdfUnique = 0;
        cdfsUnique = 1;
        pVals(1:length(sortedInt)) = 0;
    end
%     [~, cdfsEnd] = pdf_cdf_emccd(min(U,max(intUnique)+1),lambda,gain, adFactor, countOffset, roNoise,L,U);
%     cdfsUnique = cdfsUnique./cdfsEnd;
%     pdfUnique = pdfUnique./cdfsEnd;


% 
%     cdfsUnique(intUnique>floor(U)) = nan; % should not be the case
%     pdfUnique(intUnique>floor(U)) = nan;

%     [~ , cdfsUnique] = pdf_cdf_emccd(intUnique,lambda,chipPars,N);



end
