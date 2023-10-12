I = imread(filenames{1}{1}{1});


figure
tiledlayout(1,3)
nexttile
imagesc(I)
nexttile
plot(stats{1}.histAll(1:end-1));
hold on
plot(stats{1}.pdf*stats{1}.stats.nBg);
nexttile
plot(stats{1}.pval)

%%
gofVal = fzero(@(x) 0.01-chi2pdf(x,3 ),10); % 0.01 gof value
figure
tiledlayout(3,3,'TileSpacing','compact','Padding','tight')
for idx1=1:3:7
    I = imread(filenames{1}{1}{idx1});
    sampIm = mat2gray(I);
    J = imadjust(sampIm,[min(sampIm(:)) min(1,4*median(sampIm(:)))]);

    nexttile
    imshow(J,'InitialMagnification','fit');
    colormap gray
    if idx1==1
        title('a)','Interpreter','latex')
    end

    
    nexttile
    hold on
    if idx1==1
        title('b)','Interpreter','latex')
    end

    structRes = stats{idx1};
    intThreshBg =   stats{idx1}.intThreshBg;
    histAll = stats{idx1}.histAll;
    stats1 = stats{idx1}.stats;
    lambdaBg = stats{idx1}.lambdaBgMLE(intThreshBg);
%         nexttile    
    binPos = 1:structRes.LU(2) + 0.5;
    [minVal , idx] = min(abs(binPos - intThreshBg));
    h1 = bar(binPos(1:idx),histAll(1:idx),1); 
    set(h1,'FaceColor',[0.4 0.6 0.9])
    set(h1,'EdgeColor','none')
    hold on
    counts = histcounts(I,[1:max(I(:))]);
%     h2 = bar(binPos(idx+1:end),histAll(idx+1:end-1),1); 
    h2 = bar(binPos(idx+1):length(counts),counts(idx+1:end),1); 

    set(h2,'FaceColor',[1 0.5 0.3])
    set(h2,'EdgeColor','none')
    hold on
    binCountsFit = stats1.nBg.*structRes.pdf;
    % binPos = binEdges(1:end-1) + diff(binEdges)/2;
    plot(binCountsFit,'--','Color','black','LineWidth',2)
    
    % Set figure labels, etc
    if idx1==7
        xlabel('Image counts','Interpreter','latex')
        ylabel('Histogram counts','Interpreter','latex')
    end
    % set(gca,'Fontsize',15)
    % axis([30 80 0 46000])
    [a,b] = ind2sub([8 8],idx1)
%     title(['(b) Fit for tile \{',num2str(a),',',num2str(b) , '\}'],'Interpreter','latex')
    % axis equal
    pbaspect([1 0.8 0.8])
%     legendEntry = strcat(['Fit, $\lambda_{bg} =  ' num2str(lambdaBg,2) ', N_{icr}^{bg}=' num2str(intThreshBg) '$']);
%     lgnd = legend('Image counts, true background','Image counts, not true background',legendEntry,'Interpreter','latex');
%     lgnd.Layout.Tile = 'south';
    nexttile
    hold on
    if idx1==1
        title('c)','Interpreter','latex')
    end

    plot(stats{idx1}.chi2Score)
    plot([1 length(stats{idx1}.chi2Score)],[gofVal gofVal],'red')
    set(gca, 'YScale', 'log')
    if idx1==7
        xlabel('Image count truncation point','Interpreter','latex')
        ylabel('$\chi^2$ score','Interpreter','latex')
    end
end

print('output/FigS4.eps','-depsc','-r300');
