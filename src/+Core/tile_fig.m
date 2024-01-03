function [] = tile_fig(statsAll,intthreshPars,lambdaBg,matrixBlocks,idx1,aa,bb)


    nexttile
%     col
    pixelsize = 160;
 image.imAverage = matrixBlocks(:,:,idx1)
    sampIm = mat2gray(image.imAverage);
    minInt = min(sampIm(:));
    medInt = median(sampIm(:));
    maxInt = max(sampIm(:));
    J = imadjust(sampIm,[minInt min(1,4*medInt)]);
%     matrixBlocksJ = reshape(permute(reshape( double(J),T,size( J,1)/T,T,[]),[1,3,2,4]),T,T,[]);
%     out = imtile(matrixBlocksJ,'thumbnailsize',[64 64]);
%     figure,imshow(out)

%     IM3 = padarray(matrixBlocksJ,[2 2],nan,'both');
%     out = imtile(pagetranspose(IM3),'thumbnailsize',[64 64],'BorderSize', 2, 'BackgroundColor', 'cyan');
%     figure,imshow(out')

    imshow(J,'InitialMagnification','fit');
%     colormreap winter
    %imshow(images.imAverage/max(images.imAverage(:)))
    hold on    
    % scale bar (ten microns in number of pixels)
    nPixels = 1e3/pixelsize;
    x = [5, 5 + nPixels ];
    y = [0.9*size(sampIm,1) , 0.9*size(sampIm,1)];
    plot(x,y,'Linewidth',5,'Color','white')
    text(0,0.05,'1 micron','Fontsize',10,'Color','white','Units','normalized')
    title(['(',aa,')'],'Interpreter','latex')
    set(gcf, 'InvertHardCopy', 'off');


    
structRes = statsAll{idx1};
intThreshBg =   intthreshPars(idx1);
lambdaBg = lambdaBg(idx1);
histAll = statsAll{idx1}.histAll;
stats = statsAll{idx1}.stats;
 nexttile    
binPos = 1:structRes.LU(2) + 0.5;
[minVal , idx] = min(abs(binPos - intThreshBg));
h1 = bar(binPos(1:idx),histAll(1:idx),1); 
set(h1,'FaceColor',[0.4 0.6 0.9])
set(h1,'EdgeColor','none')
hold on
h2 = bar(binPos(idx+1:end),histAll(idx+1:end-1),1); 
set(h2,'FaceColor',[1 0.5 0.3])
set(h2,'EdgeColor','none')
hold on
binCountsFit = stats.nBg.*structRes.pdf;
% binPos = binEdges(1:end-1) + diff(binEdges)/2;
plot(binCountsFit,'--','Color','black','LineWidth',2)

% Set figure labels, etc
xlabel('Image counts','Interpreter','latex')
ylabel('Histogram counts','Interpreter','latex')
% set(gca,'Fontsize',15)
% axis([30 80 0 46000])
[a,b] = ind2sub([8 8],idx1)
title(['(',bb,') Fit for tile \{',num2str(a),',',num2str(b) , '\}'],'Interpreter','latex')
% axis equal
pbaspect([1 0.8 0.8])
legendEntry = strcat(['Fit, $\lambda_{bg} =  ' num2str(lambdaBg,2) ', N_{icr}^{bg}=' num2str(intThreshBg) '$']);
lgnd = legend('Image counts, true background','Image counts, not true background',legendEntry,'Interpreter','latex','Location','eastoutside');
% lgnd.Layout.Tile = 'east';
end

