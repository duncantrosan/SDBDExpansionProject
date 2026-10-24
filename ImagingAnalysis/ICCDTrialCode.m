clear 
close all
%%% Defintions
%%% Options 
LogicGif = 1; % Create GIF Output


LightImagePath = 'C:\Users\dptrosan\Documents\MATLAB\Data\2023\04_10_2023_Trial_FilimentCounting\04_10_2023_Trial_FilimentCounting\2023 April 10 15_13_22-Frame-016.csv';
path = uigetdir('C:\Users\dptrosan\Documents\MATLAB\Data\2023\04_10_2023_Trial_FilimentCounting');
FileList = dir(path);
Nfiles= length(FileList);
AllImages = zeros(1024,1024,Nfiles);
LightImage = ImageProcess(LightImagePath,1200);
j = 0;

%% Pimax settings
StartTime = 80*10^-6;
GateWidth = 80*10^-9;
TimePoints = 141;
EndTime = 160*10^-6;
t = linspace(StartTime,EndTime,TimePoints);
i = 1;


for n = 3:Nfiles 
    filename = fullfile(FileList(n).folder,FileList(n).name);
    [Image,Logic(n),FilimentCount(n),Values{n}] = ImageProcess(filename,1200);
    if Logic(n) == 1
        j = j+1;
        AllImages(:,:,j) = Image;
        TimeIndex(j) = n;
        FilimentAreaTemp = Values{1,n}{:,'A` - Measured Area [mm^2]'};
        FilimentAreaLength = length(FilimentAreaTemp);
        for k = 1:FilimentAreaLength
            AllFilimentArea(i) = FilimentAreaTemp(k);
            i = i+1;
        end
    end

end
AllImages = AllImages(:,:,1:j);
position = [400,50];
for n = 1:j

     h = figure;
    ax1 = axes;
    imagesc(ax1,LightImage);
    colormap(ax1,'gray');
    ax1.Visible = 'off'; 
    ax2 = axes;
    imagesc(ax2,AllImages(:,:,n),'alphadata',AllImages(:,:,n)>0);
    colormap(ax2,'jet');
    caxis(ax2,[min(nonzeros(AllImages(:,:,1))) max(nonzeros(AllImages(:,:,1)))]);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    colorbar;
    ImageFilename = ['C:\Users\dptrosan\Documents\MATLAB\Results\VDBDFiliment/','Image_' , num2str(n),'.png'];
    exportgraphics(h,ImageFilename)
end




