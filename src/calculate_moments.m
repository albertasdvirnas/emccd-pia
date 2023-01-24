
% Data for Fig1:
% calculate_moments("C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 beads high conc\100x\*.tif");
% calculate_moments('C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\100x\*.tif');
% calculate_moments('C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\Jason_oskar_20191125_ixon_statistics\20x\*.tif');
function name= calculate_moments(dirf)
    % function to calculate moments

    % three parameters: zoom (20x), gain (00,50,100,300) and lamp (00
    % to 100"
if nargin < 1
    files = dir("C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 beads high conc\100x\*.tif");
else
    files = dir(dirf);
end

filesC = arrayfun(@(x) fullfile(files(x).folder,files(x).name),1:length(files),'un',false);

expStruct = [];
for i=1:length(filesC)
    [fb,fd,fe] = fileparts(filesC{i});
    fieldVals = strsplit(fd,'_');
    expStruct.magnification{i} = fieldVals{1};
    expStruct.gain(i) = str2double(fieldVals{2}(5:end));
    expStruct.lamp(i) = str2double(fieldVals{3}(5:end));  
end
    

% [fb,fd,fe] = fileparts(filesC{1});

intensities = ['00',arrayfun(@(x) num2str(x), 10:10:100, 'UniformOutput', false)];

 
gain = [];
gain.means = [];
gain.vars = [];

calcBox = 0;

if calcBox
    limX = 100:400; % limit to the center,
    limY = 100:400;
end

for i = 1:length(filesC)
    A = [];
    nF = length(imfinfo(filesC{i}));

    for k=1:nF
        A(:,:,k) = imread(filesC{i},k);
    end
    
    if calcBox
        A = A(limX,limY,:);
    end 
    
%     file = filesC{cellfun(@(x) ~isempty(strfind(x,strcat(['_lamp' intensities{gk} '_']))),filesC,'UniformOutput',true)};
  
    mMat =  mean(double(A),3);
    vMat =  var(double(A),1,3);
    expStruct.means{i} = mMat(:);
    expStruct.vars{i} = vMat(:);
end

% file = filesC{cellfun(@(x) ~isempty(strfind(x,strcat(['_lamp' intensities{1} '_']))),filesC,'UniformOutput',true)};

splt =strsplit(fd,'_');

name =  strcat([splt{1} '.mat']);
save(name,'expStruct', '-v7.3');



end