function [stats,binarizedImage] = estimate_statistics(im,structRes,pValThresh)
% estimates statistics

if nargin < 3
    pValThresh = 0.01;
    im = images.imAverage;
end
% binarize image
% structRes
nPixels = numel(im);
alpha = 1;
    
% find the value where pvalue=1-cdf > pValThresh. Untruncated or truncated?
intThresh = find(1-structRes.cdf < pValThresh,1,'first'); 
if isempty(intThresh)
    intThresh = length(structRes.cdf); % take last one
end
stats.intThresh = intThresh;
% Binarize using intensity threshold
binarizedImage = imbinarize(im,intThresh);

% estimate nbg
cdfIntensities = structRes.LU:structRes.lowestIntThresh;

stats.yval =    arrayfun(@(x) sum(structRes.histAll(1:x)), cdfIntensities);
stats.xval =    arrayfun(@(x) structRes.cdf(x), cdfIntensities);
nBg =  stats.xval'\stats.yval';
stats.nBg = nBg;% (eq 12) manuscript originally

% Total number of signal pixels
nSignal = max(0,nPixels - nBg);

%  Number of signal pixels which are "black" (below the intensity
%  threshold)
nBlack = sum(structRes.histAll(1:intThresh));

% eq 13
nBlackBg = nBg*structRes.cdf(intThresh);
  
%
nBlackSignal = max(0,nBlack - nBlackBg);
nBlackSignal = min(nBlackSignal,nSignal);
 
% should be controled by pValThresh
pBlackBg = (nBlackBg + alpha)/(nBg + 2*alpha);

% ok to be 0
pBlackSignal = (nBlackSignal + alpha)/(nSignal + 2*alpha);

fBg = (nBg + alpha)/(nPixels + 2*alpha);  % ratio of black pixels
stats.misClassRate = (1-pBlackBg).*fBg + pBlackSignal*(1-fBg); % ra
    



fdr = (nBg*(1-structRes.cdf(intThresh)))/sum(structRes.histAll(intThresh+1:end));
FOR = (max(0,sum(structRes.histAll(1:intThresh))-nBg*structRes.cdf(intThresh)))/sum(structRes.histAll(1:intThresh));

stats.tnr = pBlackBg;
stats.fpr = 1 - pBlackBg;
stats.fnr = pBlackSignal;
stats.tpr = 1-pBlackSignal;
stats.fdr = fdr;
stats.for = FOR;

end

