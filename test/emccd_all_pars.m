function [chi2Score,paramEst,pdf,cdf,histAll,LU] = emccd_all_pars(I, gainGuess, adFactorGuess, roNoiseGuess)

data = sort(I(:));
m = length(data);


offsetGuess = min(I(:))+1;

lamGuess =     abs((data(round(m/2)) - offsetGuess)/(gainGuess/adFactorGuess));

% Prep data for truncated fit
import Core.calc_bounds;
[L, U, EX, STD] = calc_bounds(lamGuess, gainGuess, adFactorGuess, offsetGuess, roNoiseGuess,6);
% L = 1;
LU = [max(1,ceil(L)) floor(U)];

Nthresh = floor(U);
% Nthresh = 50;
% Get bin edges
binCenters = [max(1,ceil(L)):1:floor(U)+1 inf]; % bin edges shifted by half  

% histogram for intensities. Calculated only once!
histAll = zeros(1,binCenters(end-1));
histAll(binCenters(1:end-1)) = histcounts(data,binCenters-0.5)';

paramsGuess = [lamGuess, gainGuess, adFactorGuess, offsetGuess, roNoiseGuess];


binPosAll = max(1,ceil(L)):Nthresh;

binPos = find(histAll>1,1,'first'):Nthresh;
import Core.log_likelihood_trunc_dist_multi; 
logL = @(params,data,cens,freq,trunc) log_likelihood_trunc_dist_multi(params,data,cens,freq,trunc);

opt = statset('MaxIter',5000,'MaxFunEvals',4000,'FunValCheck','on');

% use mle with negative loglikelihood (nloglf)
[paramEst, pci] = mle(binPos,'nloglf',logL,'start',paramsGuess,'lowerBound',zeros(1,length(paramsGuess)),'Frequency',histAll(binPos),'TruncationBounds',LU,'Options',opt);

lambdaBg = paramEst(1);
gain = paramEst(2);
adFactor = paramEst(3);
offset = paramEst(4);
roNoise = paramEst(5);



cdf = nan(1,length(histAll)-1);
pdf = nan(1,length(histAll)-1);
import Core.pdf_cdf_emccd;

[pdf(binPosAll), cdf(binPosAll)] = pdf_cdf_emccd(binPosAll, lambdaBg, gain, adFactor, offset, roNoise);

% [a,b] = pdf_cdf_emccd(375, lambdaBg, gain, adFactor, offset, roNoise, LU(1), LU(2));

nQuantiles = 5;
[chi2Score] = chi2_fun(histAll, pdf,cdf,Nthresh,LU,nQuantiles);


end

