
% 1) Simulate data
% 2) Estimate lambda based on Nthresh
% 3) Calculate pdf/cdf for intensities in range
% 4) Create bins for 
dims = 8*[128 128];
numRuns = 100;
%% Chi square dist generation
%     truncationIdx = 1;      % truncation is done at truncationIdx*100 % of the data  
%     nRandNumbersBg = 1E3;     % number of random numbers generated
%     nRuns = 100;                % number of simulation runs
nQuantiles = 5;          % number of quantiles used for binning
% Nthresh = 50;

  tic
SNRVals = [3]; 
[filenames,chipPars,data] = simulate_random_beads_full(100, 20, 100, 1/800,SNRVals, 1, [],dims); % zoom 100x

gain = chipPars.gain;
adFactor = chipPars.adFactor;
roNoise = chipPars.roNoise;
offset = chipPars.countOffset;
r = gain/adFactor;

%% Estimate lambdaBg
intensities = data{1}{1}.image;
sortI = sort(data{1}{1}.image,'ascend');

intVals = min(intensities)+10:max(intensities);


lambdaBg = nan(1,length(intVals));
chi2Score = nan(1,length(intVals));
for Nthresh = intVals;
  
    try %histAll,lamGuess,Nthresh,gain,adFactor,offset,roNoise, LU, opt
        [lambdaBg(Nthresh), pci, L, U] = est_lambda(sortI, Nthresh, gain, adFactor, offset, roNoise, r);
        [chi2Score(Nthresh)] = chi2_fun(sortTruncI,Nthresh,L,U,lambdaBg(Nthresh),gain,adFactor,offset,roNoise,nQuantiles);

    catch
    end;
% lambdaBg

end

pval = chi2pdf(chi2Score,nQuantiles-2 ); 
figure,plot(pval)

% chi2Score (chi2Score < )

figure,plot(intVals,chi2Score(intVals));hold on
plot(intVals,lambdaBg(intVals))
legend({'$\chi^2$ score','$\lambda_{bg}$'},'Interpreter','latex')
title(['Nquantiles = ', num2str(nQuantiles)])