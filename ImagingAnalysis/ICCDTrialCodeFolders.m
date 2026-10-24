clear 
close all
%%% Defintions
%%% Options 
LogicGif = 1; % Create GIF Output
FindConversionFactor = 'yes';

LightImagePath = 'C:\Users\dptrosan\Documents\MATLAB\Data\2023\04_10_2023_Trial_FilimentCounting\04_10_2023_Trial_FilimentCounting\2023 April 10 15_13_22-Frame-016.csv';
path = uigetdir('C:\Users\dptrosan\Documents\MATLAB\Data\2023\04_10_2023_Trial_FilimentCounting');
FolderList = dir(path);
Nfolders= length(FolderList);
LightImage = ImageProcess(LightImagePath,0,1);
BackgroundImage = ImageProcess(LightImagePath,0,1);
j = 0;



%% Find ConversionFactor if Logic Conversion Factor is set to 'yes'
if strcmp(FindConversionFactor,'yes') == 1
    figure;
    imshow(LightImage,[1,max(max(LightImage))])
    for n = 1:5
        disp(['Please indentify two distance points ',num2str(5-(n-1)),' more times'])
    [x,y] = ginput(2);
    Pixels(n) = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    disp(['Pixel length = ',num2str(Pixels)])
    end
    Distance = input('Distance in mm = ');
    ConversionFactor = mean(Distance./Pixels); %% Conversionfactor in mm/pixel
    save('ConversionFactorICCD.mat','ConversionFactor')
elseif strcmp(FindConversionFactor,'no') == 1
    load('ConversionFactorICCD.mat')
    %ConversionFactor = ConversionFactor
else
    disp('Warning Please Put in either yes or no for the value of FindConversionFactor')
    return
end


%% Pimax settings
StartTime = 80*10^-6;
GateWidth = 80*10^-9;
TimePoints = 143;
EndTime = 160*10^-6;
t = linspace(StartTime,EndTime,TimePoints);
i = 1;
AllImages =zeros(1024,1024,143);
for n = 3:Nfolders
    Foldername = fullfile(FolderList(n).folder,FolderList(n).name);
    FileList = dir(Foldername);
    Nfiles = length(FileList);
    disp(Foldername)
    j = 1;
    for nn = 3:Nfiles
        filename = fullfile(FileList(nn).folder,FileList(nn).name);
        [Image,Logic(n,nn),FilimentCount(n,nn),Values{n,nn},Location{n,nn}] = ImageProcess(filename,700,ConversionFactor);
        if Logic(n,nn) == 1
            j = j+1;
            if n == Nfolders
            AllImages(:,:,j) = Image;
            end
            TimeIndex(n,j) = nn;
        end

    end
    disp(['Finished ',FolderList(n).name])
end


i = 1;
g = 1;
for n = 3:Nfolders
    for nn = 3:Nfiles
        if Logic(n,nn) == 1
            FilimentAreaTemp = Values{n,nn}{:,'A` - Measured Area [mm^2]'};
            EccentricityTemp = Values{n,nn}{:,'Eccentricity'};
            EccentricutyLength = length(EccentricityTemp);
            IntensityAverageTemp = Values{n,nn}{:,'Average Intensity'};
            IntensityTemp =  Values{n,nn}{:,'Intensity'};
            RadialPositionTemp =  Values{n,nn}{:,'Radial Position / mm'};
            FilimentAreaLength = length(FilimentAreaTemp);
            for k = 1:FilimentAreaLength
                AllFilimentArea(i) = FilimentAreaTemp(k);
                AllEccentricity(i) = EccentricityTemp(k);
                AllIntensity(i) =IntensityTemp(k);
                IntensityAverage(i) = IntensityAverageTemp(k);
                AllRadialPosition(i) = RadialPositionTemp(k);
                if IntensityTemp(k) >= 50000 && IntensityTemp(k) <= 1250000
                    IntensitySingle(g) = IntensityTemp(k);
                    FilimentAreaSingle(g) = FilimentAreaTemp(k);
                    EccentricitySingle(g) = EccentricityTemp(k);
                    IntensityAverageSingle(g) = IntensityAverageTemp(k);
                    g = g+1;
                end
                i = i+1;
            end
            
        end
    end
end
AverageFilimentCount = mean(FilimentCount);



AllImages = AllImages(:,:,1:j);
position = [400,50];
AllImages(AllImages<=700) =0;
for n = 2:j

     h = figure;
    ax1 = axes;
    imagesc(ax1,LightImage);
    colormap(ax1,'gray');
    ax1.Visible = 'off'; 
    ax2 = axes;
    imagesc(ax2,AllImages(:,:,n),'alphadata',AllImages(:,:,n)>0);
    colormap(ax2,'jet');
    clim(ax2,[min(nonzeros(AllImages(:,:,2))) max(nonzeros(AllImages(:,:,2)))])
    [LocationLength,~] = size(Location{Nfolders,TimeIndex(Nfolders,n)});
    for i = 1:LocationLength
        if Values{Nfolders,TimeIndex(Nfolders,n)}{i,'Intensity'} <= 1250000 
            co = 'red';
        elseif  Values{Nfolders,TimeIndex(Nfolders,n)}{i,'Intensity'} > 1250000 && Values{Nfolders,TimeIndex(Nfolders,n)}{i,'Intensity'} <= 3000000
            co = 'blue';
        elseif Values{Nfolders,TimeIndex(Nfolders,n)}{i,'Intensity'} > 3000000
            co = 'yellow';
        end
    rectangle('Position', Location{Nfolders,TimeIndex(Nfolders,n)}(i,:), 'EdgeColor', co, 'LineWidth', 2); hold on
    text(Location{Nfolders,TimeIndex(Nfolders,n)}(i,1),Location{Nfolders,TimeIndex(Nfolders,n)}(i,2)-15,['Intensity = ',num2str(Values{Nfolders,TimeIndex(Nfolders,n)}{i,'Intensity'})]...
        ,'Color','white','FontSize',6);
     text(Location{Nfolders,TimeIndex(Nfolders,n)}(i,1),Location{Nfolders,TimeIndex(Nfolders,n)}(i,2)+Location{Nfolders,TimeIndex(Nfolders,n)}(i,4)+15,['Area = ',num2str(Values{Nfolders,TimeIndex(Nfolders,n)}{i,'A` - Measured Area [mm^2]'}), ' mm^2']...
        ,'Color','white','FontSize',6)
    end
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    colorbar;
    ImageFilename = ['C:\Users\dptrosan\Documents\MATLAB\Results\VDBDFiliment/','Image_' , num2str(n),'.png'];
    exportgraphics(h,ImageFilename)
end




