clc
clear all
close all
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 22;
% Check that user has the Image Processing Toolbox installed.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
  % User does not have the toolbox installed.
  message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
  reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
  if strcmpi(reply, 'No')
    % User said No, so exit.
    return;
  end
end
%%
%UI parameters
brightness = 2;
cropFrame = [2700 50 800 1830]; %Define crop frame [top-left x, top-left y, bottom-right from-x1, bottom-right, from-y1]
flatDepth = 3;
threshold = 0.5;
%   Input values measured from image-set
con = 37.6; % conversion factor in pixels per mm
xNoz = 562; % x - position for nozzle
yNoz = 1778; % y - position for nozzle
%%
%Load/Write Directories
selpath = uigetdir
cd(selpath);
mkdir crops
mkdir corrections
mkdir filled
mkdir data
cropDir = strcat(selpath,'/crops');
corDir = strcat(selpath,'/corrections');
fillDir = strcat(selpath,'/filled');
dataDir = strcat(selpath,'/data');
files1 = dir('*.csv*'); %consolidate bubble images
flats = dir('*flat*.JPG'); %consolidate flat-corrections
totalFiles = numel(files1);
totalFlats = numel(flats); 



% Definie OPtions 


        opts = delimitedTextImportOptions("NumVariables", 5);
        opts.DataLines = [1, Inf];
        opts.Delimiter = ",";
        opts.VariableNames = ["Var1", "Var2", "VarName3", "VarName4", "VarName5"];
        opts.SelectedVariableNames = ["VarName3", "VarName4", "VarName5"];
        opts.VariableTypes = ["string", "string", "double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts = setvaropts(opts, ["Var1", "Var2"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var1", "Var2"], "EmptyFieldRule", "auto");

%%
%Crop raw images
parfor i =1:totalFiles
 %  %Load all files in directory
        Fileaddress{i,1}=strcat(selpath,'/',files1(i).name);

        % Import the data
        Table = readtable(Fileaddress{i,1}, opts);
        InputMatrix = table2array(Table);

        %% Clear temporary variables

        ImageMatrix = zeros(1024,1024);
        n = 1;
        % Restructure Information from Pimax Camera
        for i = 1:1024
            for j = 1:1024
                ImageMatrix(i,j) = InputMatrix(n,3);
                n = n+1;
            end
        end
        


        file{i} = imread(Fileaddress{i,1});
    %Crop images
        img_crop{i} = imcrop(file{i},cropFrame);
        cd(cropDir) % go to dir where you want to save updated files
        writeFileName = strcat('crop_',num2str(i),'.PNG');
        imwrite(img_crop{i} ,writeFileName)
           
end
        clear opts
for i =1:totalFlats
   %Load all files in directory
       Flataddress{i,1}=strcat(selpath,'/',flats(i).name);
       flat{i} = imread(Flataddress{i,1});
    %Crop images
        flat_crop{i} = imcrop(flat{i},cropFrame);
    %Write Cropped images
            cd(cropDir) % go to dir where you want to save updated files
            writeFileName = strcat('flat_',num2str(i),'.PNG');
            imwrite(flat_crop{i} ,writeFileName)
            
end
%%
%Image Processing
parfor i = 1:totalFlats
    flatMat{i} = im2double(flat_crop{i});
end

for i = 1:totalFiles
   imgMat{i} = im2double(img_crop{i});
   numerator = 256 .* imgMat{i};
        for j = 1:flatDepth
        denominator = flatMat{j} + 1;
        cor = numerator ./ denominator;
        numerator = 256 .* cor;
        end
    corHold{i} = mat2gray(cor);
            cd(corDir) % go to dir where you want to save updated files
            writeFileName = strcat('correction_',num2str(i),'.PNG');
            imwrite(corHold{i} ,writeFileName)
    corHold{i} = imread(writeFileName);
    correction{i} = immultiply(corHold{i},brightness);
    imwrite(correction{i} ,writeFileName)
end

%%
parfor i = 1:totalFiles
grayImage{i} = rgb2gray(correction{i});
invertImage{i} = imcomplement(grayImage{i});
binaryImage{i} = imbinarize(invertImage{i},threshold);
binaryImage{i} = bwareaopen(binaryImage{i}, 100);
binaryImage{i} = imfill(binaryImage{i}, 'holes');
[labeledImage{i}, bubbleCount{i}] = bwlabel(binaryImage{i});
colorLabeledImage{i} = label2rgb (labeledImage{i}, 'hsv', 'k', 'shuffle');
            cd(fillDir) % go to dir where you want to save updated files
            writeFileName = strcat('filled_',num2str(i),'.PNG');
            imwrite(colorLabeledImage{i} ,writeFileName)
end
%%
% Measure Images
parfor i = 1:totalFiles
stats{i} = regionprops('table',binaryImage{i},'Area','Centroid',...
    'MajorAxisLength','MinorAxisLength',...
    'Eccentricity','Orientation');

tranPos{i} = stats{i}{:,{'Centroid'}}; %load bubble position
diffPos{i} = [xNoz.*ones(size(tranPos{i},1),1),yNoz.*ones(size(tranPos{i},1),1)]; % create matrix of nozzle positions per bubble in image
tranPos{i} = abs(tranPos{i}-diffPos{i}); % transform reference frame to nozzle
convert{i} = [tranPos{i}/con,stats{i}{:,{'MajorAxisLength','MinorAxisLength'}}/2/con,stats{i}{:,{'Area'}}/con^2]; % convert bubble values from pixels to cm and sq.cm
calcArea{i} = convert{i}(:,3).*convert{i}(:,4)*pi; % calculate area from axes
calcVolume{i} = 4.*calcArea{i}.*convert{i}(:,3); % calculate volume ASSUMING bubbles are oblate ellipsoid
ellipse{i} = abs(100.*minus(calcArea{i},convert{i}(:,5))./calcArea{i}); % percent difference between bubble and absolute ellipse
stats2{i} = [convert{i} calcArea{i} calcVolume{i} ellipse{i} stats{i}{:,{'Eccentricity','Orientation'}}];% new table of calculated values per bubble
digits(4);
values{i} = array2table(stats2{i});
values{i}.Properties.VariableNames = {'x - Position [mm]' 'y - Position [mm]'...
    'a - Major Axis [mm]' 'b - Minor Axis [mm]' 'A` - Measured Area [mm^2]'...
    'A - Calc Area [mm^2]' 'V - Calc Volume [mm^3]' 'Ellipse Deviation [%]'...
    'Eccentricity' 'α - Angle [ ˚ ]'}
    cd(dataDir)
    writeFileName = strcat('data_',num2str(i),'.csv');
writetable(values{i},writeFileName)
end

