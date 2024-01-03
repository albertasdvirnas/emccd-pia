function [outputArg1,outputArg2] = figS11(inputArg1,inputArg2)
% https://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/

%%
I = imread('https://blogs.mathworks.com/images/steve/60/nuclei.png');
I_cropped = I(400:400+511, 465:465+511);
% imshow(I_cropped)
% imcredit('Image courtesy of Dr. Ramiro Massol')


I_eq = adapthisteq(I_cropped);
imshow(I_eq)

bw = im2bw(I_eq, graythresh(I_eq));
imshow(bw)

bw2 = imfill(bw,'holes');
bw3 = imopen(bw2, ones(5,5));
bw4 = bwareaopen(bw3, 40);
bw4_perim = bwperim(bw4);
overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
imshow(overlay1)

I = I_eq;

I(~bw4) = nan;

%%
rng('default')

SNRVals =3:0.5:10;  

im = double(I)./double(max(I(:)));
im(im~=0) = im(im~=0)+1;

filenames = simulate_image_full(100, 20, 100, 1,SNRVals, 1, 'synthSingleFrameImage',[512 512], im,[]); % zoom 100x



outFig3 = 'output\Fig2S.eps';
chipParsCur = Core.chippars_specific_gain(chipParsSAll,3);


chipParsCur.pval = 0.01;
% emccdpia_thresholding(filenames{1}{1}(6),SNRVals,chipParsCur,outFig3,chipParsCur.pval,1);

[results,stats] = emccdpia_thresholding(filenames{1}{1}(1:end),SNRVals,chipParsCur,outFig3,chipParsCur.pval,1);


figure_sup


figS10_segmentation(chipParsSAll, 'output\Fig2Sc.eps')

%(imageFilenames,

%%

end

