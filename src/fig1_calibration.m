function chipPars = fig1_calibration(no_gain_mat, gain_mat,outFig)
    
    %   Args:
    %       no_gain_mat : file with no gain statistics
    %       gain_mat : mat file with gain statitistics 
    %
    %   Returns:
    %       chipPars : estimated chip parameters
    if nargin < 1
        no_gain_mat = '100x.mat';
        gain_mat = '20x.mat';
    end

    %
    dataGain0 =  importdata(no_gain_mat); % gives the same
    dataGain =  importdata(gain_mat); % why 100 for offset estimation?


    idx = find(dataGain0.gain==0); % indices for gain
    int = dataGain0.lamp(idx); % intensities
    [sortvals,sortidx] = sort(int,'asc');

    gain = cell(1,4);
    gain{1}.means = cell2mat(dataGain0.means(idx(sortidx)));
    gain{1}.vars = cell2mat(dataGain0.vars(idx(sortidx)));
    
    % remove first / this could be used separately for offset/sigma
    % estimation
    gain{1}.means = gain{1}.means(2:end,:);
    gain{1}.vars = gain{1}.vars(2:end,:);

    gval = [50 100 300];
    for ii=1:length(gval)
        idx = find(dataGain.gain == gval(ii)); % indices for gain
        intv = dataGain.lamp(idx); % intensities
        %
        [sortvals,sortidx] = sort(intv,'asc');
        gain{ii+1}.means = cell2mat(dataGain.means(idx(sortidx)))';
        
        gain{ii+1}.means = gain{ii+1}.means(2:end,:);

        gain{ii+1}.means = gain{ii+1}.means(:);

        gain{ii+1}.vars = cell2mat(dataGain.vars(idx(sortidx)))';
        
        gain{ii+1}.vars = gain{ii+1}.vars(2:end,:);

        gain{ii+1}.vars  = gain{ii+1}.vars(:);
        
    end

%
% plot mean/var should be increasing
% ix=1
%  figure,plot(gain{ix}.means(:),gain{ix}.vars(:),'x')
% figure,plot(arrayfun(@(x) mean(gain{ix}.means(x,:)),1:size(gain{ix}.means,1)),arrayfun(@(x) mean(gain{1}.vars(x,:)),1:size(gain{ix}.means,1)))

% gain = importdata('meanvardatasim.mat');

% remove 0.1% of  smallest and 0.1% largest values
% cuts = [0.01 .99];
cuts = [0.01 0.99];


varcalc = cell(1,length(gain));
mred = cell(1,length(gain));
chipPars= [];
chipPars.inid = [];
N = 1000; % over how many points
Ntrials = 2500; % number of trials (for estimating confidence) /max Ntrials

for i=1:length(gain)
    
    % first calculate analog to digital factor f
    mVec = gain{i}.means(:); % means
    varVec = gain{i}.vars(:);
    
    idx = ceil(numel(mVec) * cuts);
    [sortM,sIdx] = sort(mVec);
    sortV = varVec(sIdx);
    mVecRed = sortM(max(idx(1),1):idx(2));
    varVecRed = sortV(max(idx(1),1):idx(2));
    
    if Ntrials*N>length(mVecRed)
        warning('too many trials, not enough data points');
        Ntrials =  floor(length(mVecRed)/N);
    end
    
    pos   = randperm(length(mVecRed));
    for j=1:Ntrials
%               
        mVecTemp = mVecRed(pos((j-1)*(N)+1:j*N));
        varVecTemp = varVecRed(pos((j-1)*(N)+1:j*N));
        [mVecTemp,pos2] = sort(mVecTemp);
        varVecTemp = varVecTemp(pos2);

        [chipPars,varcalc,mred ] = calc_chip_params(i,j,varVecTemp,mVecTemp,chipPars,varcalc,mred ); 
    end
    
    if i==1 % take means of offset and adfactor for more accurate statistics
        chipPars.adFactor = 1/mean(chipPars.slope);
        chipPars.offsetFactor = mean(chipPars.offset);
    else
        % calculate count offset and readout noise
        offset_fun = @(f,g,b1,b2) f*(b1-b2)/(2*g-1);
        chipPars.countOffsetAvg(i) = offset_fun(chipPars.adFactor,mean(chipPars.gain{i}),chipPars.offsetFactor,mean(cellfun(@(x) x(1),chipPars.gaincoeffs{i})));
        % function for sigma from mean/variance relations
        sigma_fun = @(b1,delta,f) sqrt(b1-1/12+delta/f);
        
        chipPars.roNoiseAvg(i) = sigma_fun(chipPars.offsetFactor,chipPars.countOffsetAvg(i),chipPars.adFactor);
    end

        
    chipPars.inid.mVecTemp{i} =  mVecRed(pos(1:N));
    [chipPars.inid.mVecTemp{i} , pos2] = sort(chipPars.inid.mVecTemp{i} );
    chipPars.inid.varVecTemp{i} =  varVecRed(pos(1:N));
    chipPars.inid.varVecTemp{i}  = chipPars.inid.varVecTemp{i} (pos2);

end

% Figure 1: calibration: plot histograms for gain parameter estimation
plot_line_fits_calibration(chipPars,varcalc)
plot_histograms_calibration(chipPars,varcalc)
print(outFig,'-depsc','-r300')


% chipPars to latex table
offset_fun = @(f,g,b1,b2) f.*(b1-b2)./(2*g-1);
% function for sigma from mean/variance relations
sigma_fun = @(b1,delta,f) sqrt(b1-1/12+delta./f);

% chipPars.slope is 1/f
% chipPars.offset b1
% b2 is  gainPars = cellfun(@(x) x(1),chipPars.gaincoeffs{2});
% g is  (gainPars./(chipPars.slope(perm)))/2;

deltaQ = zeros(4,3);
sigmaQ = zeros(4,3);
gainQ = zeros(4,3);
for ii=2:4
% now we get delta: 50
    gainPars = cellfun(@(x) x(2),chipPars.gaincoeffs{ii});
    Cgain = cellfun(@(x) x(1),chipPars.gaincoeffs{ii});
    perm = randperm(length(gainPars));
    g = (gainPars./(chipPars.slope(perm(1:length(gainPars)))))/2;

    offsetVals = offset_fun(1./chipPars.slope((1:length(gainPars))),g,chipPars.offset(perm(1:length(gainPars))),Cgain);

    deltaQ(ii,:) = quantile(offsetVals,3);
    sigmaVals = sigma_fun(chipPars.offset(perm(1:length(gainPars))),offsetVals,1./chipPars.slope(perm(1:length(gainPars))));
    sigmaQ(ii,:) = quantile(sigmaVals,3);
    gainQ(ii,:) = quantile(g,3);

end

aduQ = quantile(1./chipPars.slope,3);

chipPars.gainQ = gainQ;
chipPars.sigmaQ = sigmaQ;
chipPars.aduQ = aduQ;
chipPars.deltaQ = deltaQ;

% 
% % fprintf('Offset with no gain: %.2f.\n',noGainCoeffs(1));
% figure(1)
% % plotpoints = randperm(numel(mVecRed),1e4);
% scatter(  chipPars.inid.mVecTemp{1}, chipPars.inid.varVecTemp{1})
% hold on
% % varsCalc = mPad * noGainCoeffs;
% plot( chipPars.inid.mVecTemp{1},varcalc{1}{1},'--','LineWidth',1,'Color','black')
% hold off
% %xlim([0 500])
% %ylim([0 30])
% xlabel('Mean image count','Interpreter','latex','FontSize',15)
% ylabel('Variance','Interpreter','latex','FontSize',15)
% fig=gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 7 5];
% set(gca,'Fontsize',15)

%chipPars.gain = gainCoeffs(2)*chipPars.adFactor/2;
fprintf('Estimated gain at  50: %.3f, std %.3f\n',mean(chipPars.gain{2}),std(chipPars.gain{2}));
fprintf('Estimated gain at 100: %.3f, std %.3f.\n',mean(chipPars.gain{3}),std(chipPars.gain{3}));
fprintf('Estimated gain at 300: %.3f, std %.3f.\n',mean(chipPars.gain{4}),std(chipPars.gain{4}));
fprintf('Estimated count offset at  50: %.2f, std %.3f.\n',mean(chipPars.countOffset{2}),std(chipPars.countOffset{2}));
fprintf('Estimated count offset at 100: %.2f, std %.3f.\n',mean(chipPars.countOffset{3}),std(chipPars.countOffset{3}));
fprintf('Estimated count offset at 300: %.2f, std %.3f.\n',mean(chipPars.countOffset{4}),std(chipPars.countOffset{4}));
fprintf('Estimated readout noise at  50: %.2f, std %.3f.\n',mean(chipPars.roNoise{2}),std(chipPars.roNoise{2}));
fprintf('Estimated readout noise at 100: %.2f, std %.3f.\n',mean(chipPars.roNoise{3}),std(chipPars.roNoise{3}));
fprintf('Estimated readout noise at 300: %.2f, std %.3f.\n',mean(chipPars.roNoise{4}),std(chipPars.roNoise{4}));

end

function [chipPars,varcalc,mred ] = calc_chip_params(i,j,varVecRed,mVecRed,chipPars,varcalc,mred  )
    %
    % function for offset from the mean/variance equations
    offset_fun = @(f,g,b1,b2) f*(b1-b2)/(2*g-1);

    % function for sigma from mean/variance relations
    sigma_fun = @(b1,delta,f) sqrt(b1-1/12+delta/f);
    
    % calculates chip params for each case
    mPad =[ones(length(mVecRed),1) mVecRed]; % will solve AX+Y = B

    gaincoeffs = mPad\varVecRed; % X = A\B  is solution to A*X = B. So we solve Y=B by first row and AX=B by second row 
%     gaincoeffs(1)= max(gaincoeffs(1),0);
%     polyfit(mVecRed,varVecRed,1)

    varcalc{i}{j} = mPad * gaincoeffs;
    mred{i}{j} = mVecRed;
    chipPars.gaincoeffs{i}{j} = gaincoeffs; 
    if i == 1
        chipPars.offset(j) = gaincoeffs(1);
        chipPars.slope(j) = gaincoeffs(2); % ADU factor f is 1/slope

%         chipPars.adFactor(j) = 1/gaincoeffs(2); % ADU factor f is 1/slope
%         fprintf('Estimated slope: %.2f.\n', chipPars.slope(j));
    else % adFactor and  offsetFactor are mean(offset) and 1/mean(slope)
        % current gain is the slope times f/2
        chipPars.gain{i}(j) = gaincoeffs(2)*chipPars.adFactor/2; 
        chipPars.countOffset{i}(j) = offset_fun(chipPars.adFactor,chipPars.gain{i}(j),chipPars.offsetFactor,gaincoeffs(1));

        %         (chipPars.offset-gaincoeffs(1))*chipPars.adFactor/(chipPars.gain{i}*2-1); % should be gain/2??
        %fprintf('Estimated count offset: %.2f.\n',chipPars.countOffset100);
        chipPars.roNoise{i}(j) = sigma_fun(chipPars.offsetFactor,chipPars.countOffset{i}(j),chipPars.adFactor);
        %         sqrt((chipPars.offset+chipPars.countOffset{i}/chipPars.adFactor)*chipPars.adFactor^2);
    end
end

function plot_line_fits_calibration(chipPars,varcalc)

figure
tiledlayout(2,2,'TileSpacing','tight')
ax1 = nexttile
% clf
% scatter(  chipPars.inid.mVecTemp{1}, chipPars.inid.varVecTemp{1})
hold on

p2 = scatter(  chipPars.inid.mVecTemp{2}, chipPars.inid.varVecTemp{2})
p3 = scatter(  chipPars.inid.mVecTemp{3}, chipPars.inid.varVecTemp{3})
p4 = scatter(  chipPars.inid.mVecTemp{4}, chipPars.inid.varVecTemp{4})

% plot( chipPars.inid.mVecTemp{1},varcalc{1}{1},'--','LineWidth',1,'Color','black')
plot( chipPars.inid.mVecTemp{2},varcalc{2}{1},'--','LineWidth',1,'Color','black')
plot( chipPars.inid.mVecTemp{3},varcalc{3}{1},'--','LineWidth',1,'Color','black')
plot( chipPars.inid.mVecTemp{4},varcalc{4}{1},'--','LineWidth',1,'Color','black')


% print gain 50 and 300 values too
hold off
xlabel('Mean image count','Interpreter','latex')
ylabel('Variance','Interpreter','latex')
title('a) Gain 50,100,300','Interpreter','latex')
% set(gca,'Fontsize',15)
%ylim([0 400])
%xlim([0 160])
% axes('Position',[0.2 0.6 0.25 0.25])
ax2 = nexttile
box on
p1 = scatter(  chipPars.inid.mVecTemp{1}, chipPars.inid.varVecTemp{1},'black')
hold on
p5 = plot( chipPars.inid.mVecTemp{1},varcalc{1}{1},'--','LineWidth',1,'Color','black')
hold off
xlabel('Mean image count','Interpreter','latex')
%xlim([0 500])
%ylim([0 30])
%xlabel('Mean image count','Interpreter','latex','FontSize',12)
%ylabel('Variance','Interpreter','latex','FontSize',12)
title('b) Gain 0','Interpreter','latex')
lgnd = legend([p1 p2 p3 p4 p5], {'Gain 0', 'Gain 50','Gain 100','Gain 300','Linear fits'},'Interpreter','latex');
set(lgnd,'color','none');
lgnd.Layout.Tile = 'east';

% fig=gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 5 3.6];


end

function fig=plot_histograms_calibration(chipPars,varcalc)

% figure(3)
% % tiledlayout(2,2,'Padding','tight')
nexttile
% clf
% scatter(  chipPars.inid.mVecTemp{1}, chipPars.inid.varVecTemp{1})
hold on

[counts1,pos] = histcounts(chipPars.gain{2},'Normalization','count');
plot((pos(2:end)+pos(1:end-1))/2,counts1,'-')
[counts2,pos] = histcounts(chipPars.gain{3},'Normalization','count');
plot((pos(2:end)+pos(1:end-1))/2,counts2,'-')
[counts3,pos] = histcounts(chipPars.gain{4},'Normalization','count');
plot((pos(2:end)+pos(1:end-1))/2,counts3,'-')
ylim([0 max([counts1 counts2 counts3])])

% legend({'Gain 50','Gain 100','Gain 300'},'Interpreter','latex')
% print gain 50 and 300 values too
hold off
xlabel('Estimated gain','Interpreter','latex')
ylabel('Histogram counts','Interpreter','latex')
title('c)','Interpreter','latex')
% set(gca,'Fontsize',15)
% xlim([0 160])
% axes('Position',[0.2 0.6 0.25 0.25])
nexttile
box on
[counts,pos] = histcounts(chipPars.slope,'Normalization','count');
plot((pos(2:end)+pos(1:end-1))/2,counts,'black-')
xlabel('$1/{\rm f}$','Interpreter','latex')

hold off
title('d) ','Interpreter','latex')

%%
% nexttile
% % clf
% % scatter(  chipPars.inid.mVecTemp{1}, chipPars.inid.varVecTemp{1})
% hold on
% 
% [counts,pos] = histcounts(chipPars.countOffset{2},'Normalization','pdf');
% plot((pos(2:end)+pos(1:end-1))/2,counts,'-')
% [counts,pos] = histcounts(chipPars.countOffset{3},'Normalization','pdf');
% plot((pos(2:end)+pos(1:end-1))/2,counts,'-')
% [counts,pos] = histcounts(chipPars.countOffset{4},'Normalization','pdf');
% plot((pos(2:end)+pos(1:end-1))/2,counts,'-')
% 
% legend({'Offset 50','Offset 100','Offset 300'},'Interpreter','latex')
% % print gain 50 and 300 values too
% hold off
% xlabel('Estimated gain','Interpreter','latex','Fontsize',15)
% ylabel('Variance','Interpreter','latex','Fontsize',15)
% set(gca,'Fontsize',15)
% ylim([0 1])
% xlim([0 160])
% axes('Position',[0.2 0.6 0.25 0.25])
% box on
% [counts,pos] = histcounts(chipPars.offset,'Normalization','pdf');
% plot((pos(2:end)+pos(1:end-1))/2,counts,'-')
% title('Offset','Interpreter','latex')

%xlim([0 500])
%ylim([0 30])
%xlabel('Mean image count','Interpreter','latex','FontSize',12)
%ylabel('Variance','Interpreter','latex','FontSize',12)
fig=gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 3.6];


end
