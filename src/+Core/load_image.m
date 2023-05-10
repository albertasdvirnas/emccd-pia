function [images] = load_image(filename)

    % Images and associated information                             
    im = imread(filename); 
    images.imAverage = double(im);
    images.registeredIm{1} = double(im);
    images.imageName = 'BeadsOnSurface';
    images.imageNumber = 1;
    images.runNo = 1;

end

