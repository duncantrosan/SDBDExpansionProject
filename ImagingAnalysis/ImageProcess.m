function [Image,Logic,FilimentNumber,Values,Location] = ImageProcess(filename,Threshold,ConversionFactor,BackgroundImage,Cropping)
% Desciption

%% Examples 

%% Readin Warnings 


XCenter = 619;
YCenter = 505;
Thresh = 0.08;
con = ConversionFactor^-1;
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

% Import the data
Table = readtable(filename, opts);
InputMatrix = table2array(Table);

%% Clear temporary variables
clear opts
ImageMatrix = zeros(1024,1024);
n = 1;
% Restructure Information from Pimax Camera
for i = 1:1024
    for j = 1:1024
        ImageMatrix(i,j) = InputMatrix(n,3);
        n = n+1;
    end
end


if Threshold == 0 
    Image = ImageMatrix;
    return
end

% Is there a filiment
% Compares the absoulte max of the file with the threshold if there is no
% filiment breaks out of function

Logic = 1;
FilimentNumber = 1;

%% Image Preprocessing Cropping,Rotating etc.
% ind2rbg 
% scaleing = turbo
% Cropping

% Background Subtraction
if exist('BackgroundImage') == 1
ImageMatrix = ImageMatrix - BackgroundImage;
end



if exist('Cropping') == 1
    
end


if max(max(ImageMatrix)) <= Threshold
    Logic = 0;
    Image = 0;
    FilimentNumber = 0;
    Values = 0;
    Location = 0;
    return
end




%% Bubble Analysis code 
OriginalMatrix = ImageMatrix;
ImageMatrix(ImageMatrix < Threshold) = 0;
BinaryImage = imbinarize(ImageMatrix,Thresh);
BinaryImage = bwareaopen(BinaryImage, 100);
BinaryImage = imfill(BinaryImage, 'holes');
[LabeledImage, FilimentNumber] = bwlabel(BinaryImage);
ColorLabeledImage = label2rgb (LabeledImage, 'hsv', 'k', 'shuffle');

            % cd(fillDir) % go to dir where you want to save updated files
            % writeFileName = strcat('filled_',num2str(i),'.PNG');
            % imwrite(ColorLabeledImage{i} ,writeFileName)


if FilimentNumber >= 1

Stats = regionprops('table',BinaryImage,'Area','Centroid',...
    'MajorAxisLength','MinorAxisLength',...
    'Eccentricity','Orientation','BoundingBox');

TranPos = Stats{:,{'Centroid'}}; %load bubble position
%diffPos{i} = [xNoz.*ones(size(tranPos{i},1),1),yNoz.*ones(size(tranPos{i},1),1)]; % create matrix of nozzle positions per bubble in image
%tranPos{i} = abs(tranPos{i}-diffPos{i}); % transform reference frame to nozzle
Convert = [TranPos/con,Stats{:,{'MajorAxisLength','MinorAxisLength'}}/2/con,Stats{:,{'Area'}}/con^2]; % convert bubble values from pixels to cm and sq.cm
CalcArea = Convert(:,3).*Convert(:,4)*pi; % calculate area from axes
CalcVolume = 4.*CalcArea.*Convert(:,3); % calculate volume ASSUMING bubbles are oblate ellipsoid
Ellipse = abs(100.*minus(CalcArea,Convert(:,5))./CalcArea); % percent difference between bubble and absolute ellipse
% Integrate over originally area 
% Find the Intesity within the area of a given region
Rect = round(Stats{:,{'BoundingBox'}});
Mask = false(size(ImageMatrix));
for n = 1:height(Stats)
    Intensity(n,1) = sum(sum(OriginalMatrix(Rect(n,1):Rect(n,1)+Rect(n,3),Rect(n,2):Rect(n,2)+Rect(n,4))));
    AverageIntensity(n,1) = Intensity(n,1)/(Rect(n,3)*Rect(n,4));
    RadialPosition(n,1) = sqrt((TranPos(n,1) - XCenter)^2 + (TranPos(n,2) - YCenter)^2 )/con ;

end

StatsLength = [Convert CalcArea CalcVolume Ellipse Stats{:,{'Eccentricity','Orientation'}} Intensity, AverageIntensity,RadialPosition];% new table of calculated values per bubble
% Find the Intesity within the area of the  given 


Values = array2table(StatsLength);
Values.Properties.VariableNames = {'x - Position [mm]' 'y - Position [mm]'...
    'a - Major Axis [mm]' 'b - Minor Axis [mm]' 'A` - Measured Area [mm^2]'...
    'A - Calc Area [mm^2]' 'V - Calc Volume [mm^3]' 'Ellipse Deviation [%]'...
    'Eccentricity' 'α - Angle [ ˚ ]','Intensity' , 'Average Intensity','Radial Position / mm'};

Location = Rect;

elseif FilimentNumber == 0;
    Values = 0;
    Logic = 0;
    Location = 0;
end


    %cd(dataDir)
    %writeFileName = strcat('data_',num2str(i),'.csv');
%writetable(values,writeFileName)



Image = ImageMatrix;
