clear
close all
%%% Defintions
%%% Options
LogicGif = 1; % Create GIF Output
FindConversionFactor = 'no';
%% Set up how data is saved
Nreplicate = 1 ;% Averages 3 files into one Ozone density picture
TimeIndex = [1000,2000,3000]; % Hz
Nset = 103;
L = 50; % mm
LogicReguclarPlotting = 0;
LogicGif = 1;
FilterFactor391 = 0.4134;
FilterFactor337 = 0.2815;
Factor391Percentage = 0.4533;
%Factor391Percentage = 0.2;
NLightImages = 3; % Number of Images in each lightImage file
%% Load in Filter Transmitance Data
TransmittanceData = load('FilterTransmitance.mat');
% Column 1 = Wavelength
% Column 2 = 253 Transmitance Column 3 = 253 error
% Column 4 = 337 Transmitance Column 5 = 337 error
% Column 6 = 380 Transmitance Column 7 = 380 error
% Column 8 = 390 Transmitance Column 3 = 390 error
% Column 2 = 700 Transmitance Column 3 = 700 error

% 
Transmittance337 = table2array(TransmittanceData.FilterTransmitance(2:end,4));
Transmittance391 = table2array(TransmittanceData.FilterTransmitance(2:end,8));
TransmittanceWavelength = table2array(TransmittanceData.FilterTransmitance(2:end,1));

figure; 
plot(TransmittanceWavelength,Transmittance337); hold on
plot(TransmittanceWavelength,Transmittance391)


Threshold = 1;
ThresholdColumns = 14;
AbelTransformLine = 14.2;
LightImagePath = 'C:\Users\dptrosan\Documents\MATLAB\Data\2023\Duncan_10_6-20231006T211740Z-001\Duncan_10_6\2023 October 06 16_47_36-Frame-1';
path337 = uigetdir('C:\Users\dptrosan\Documents\MATLAB\Data\2023\FullImagingDataSet\Fulldata\SingleAccumuationImaging','Please Select 337 Filter Directory');
path391 = uigetdir('C:\Users\dptrosan\Documents\MATLAB\Data\2023\FullImagingDataSet\Fulldata\SingleAccumuationImaging','Please Select 391 Filter Directory');
FileList337 = dir(path337);
FileList391 = dir(path391);
Nfiles337= length(FileList337);
Nfiles391 = length(FileList391);
j = 0;
BackgroundFilename = 'C:\Users\dptrosan\Documents\MATLAB\Data\2023\FullImagingDataSet\Fulldata\12_19_23Meanenergy\SmallHexKapton\FullPeriod\Background\337filterFullPeriodBG_realBG 2023 December 19 13_16_41';
BackgroundFilename391 = 'C:\Users\dptrosan\Documents\MATLAB\Data\2023\FullImagingDataSet\Fulldata\12_19_23Meanenergy\SmallHexKapton\FullPeriod\Background\391filterFullPeriodBG 2023 December 19 13_22_09';

Sigma = 11.5 *10^-18; % cm^2
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
TotalFiles = Nfiles337+Nfiles391;
PercentDone = 0;
% Readin Sample VC curve
filenameVC = 'C:\Users\dptrosan\Documents\MATLAB\Data\2023\Duncan_ICCDwithFilters\kapton337full_000_ALL';

opts1 = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts1.DataLines = [12, Inf];
opts1.Delimiter = ",";

% Specify column names and types
opts1.VariableNames = ["TIME", "CH1", "CH2"];
opts1.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts1.ExtraColumnsRule = "ignore";
opts1.EmptyLineRule = "read";

% Import the data
kapton337full000ALL = readtable(filenameVC, opts1);

%% Convert to output type
VCData = table2array(kapton337full000ALL);

%% Clear temporary variables
clear opts1

% findBackgroundImage
Table = readtable(BackgroundFilename, opts);
InputMatrix = table2array(Table);
Background = zeros(1024,1024);
n = 1;
% Restructure Information from Pimax Camera
for i = 1:1024
    for j = 1:1024
        if InputMatrix(n,1) >= Threshold
            Background(i,j) = InputMatrix(n,1);
        end
        n = n+1;
    end
end

% findBackgroundImage for 391 line
Table = readtable(BackgroundFilename391, opts);
InputMatrix = table2array(Table);
Background391 = zeros(1024,1024);
n = 1;
% Restructure Information from Pimax Camera
for i = 1:1024
    for j = 1:1024
        if InputMatrix(n,1) >= Threshold
            Background391(i,j) = InputMatrix(n,1);
        end
        n = n+1;
    end
end

for nn = 3:Nfiles337

    %% Clear temporary variables
    % Update Light Image for each data set

    filename = fullfile(FileList337(nn).folder,FileList337(nn).name);
    Table = readtable(filename, opts);
    InputMatrix = table2array(Table);
    ImageMatrix = zeros(1024,1024);
    n = 1;
    % Restructure Information from Pimax Camera
    for i = 1:1024
        for j = 1:1024
            if InputMatrix(n,1) >= Threshold
                ImageMatrix(i,j) = InputMatrix(n,1);
            end
            n = n+1;
        end
    end
    Data337All(nn-2,:,:) = (ImageMatrix - Background)/FilterFactor337;
    PercentDone = PercentDone +1;
    disp(['Percent Processed ',num2str(round(100*PercentDone/TotalFiles)),'%'])
end

for nn = 3:Nfiles391

    %% Clear temporary variables
    % Update Light Image for each data set

    filename = fullfile(FileList391(nn).folder,FileList391(nn).name);
    Table = readtable(filename, opts);
    InputMatrix = table2array(Table);
    ImageMatrix = zeros(1024,1024);
    n = 1;
    % Restructure Information from Pimax Camera
    if nn >= 18
        rows = 1000;
    else
        rows = 1024;
    end
    for i = 1:rows
        for j = 1:1024
            if InputMatrix(n,1) >= Threshold
                ImageMatrix(i,j) = InputMatrix(n,1);
            end
            n = n+1;
        end
    end
    Data391All(nn-2,:,:) = (ImageMatrix - Background391)/FilterFactor391*Factor391Percentage;
    PercentDone = PercentDone +1;
    disp(['Percent Processed ',num2str(round(100*PercentDone/TotalFiles)),'%'])
end

% OfferEoN = 244.266+7053.353*Ir-23732.444*Ir.^2;

ReducedElectricFieldData = GetData;
x = ReducedElectricFieldData(:,6);
y =  ReducedElectricFieldData(:,1);
Coeff = polyfit(x,y,4);
Fit = x.^4.*Coeff(1) + x.^3.*Coeff(2) + x.^2.*Coeff(3) + x.*Coeff(4) + Coeff(5);


for nn = 1:Nfiles391-2
    k = 0;
    EoNData = zeros(1024,1024);
    A = squeeze(Data337All(nn,:,:));
    B = squeeze(Data391All(nn,:,:));
    Ratio = B./A;
    Ratio(Ratio>= 0.6) = 0;
    % Ratio = medfilt2(Ratio,[3,3]);
    for i = 1:1024
        for j = 1:1024
                    if B(i,j) >= 800
                    % EoNData(i,j) = 244.266+7053.353*Ratio(i,j)-23732.444*Ratio(i,j)^2;
                    EoNData(i,j) = Ratio(i,j).^4.*Coeff(1) + Ratio(i,j).^3.*Coeff(2) + Ratio(i,j).^2.*Coeff(3) + Ratio(i,j).*Coeff(4) + Coeff(5);
                    k = k+1;
                    end

        end
    end
    EoNDataAll(nn,:,:) = EoNData;
    figure;
    image(B,"CDataMapping","scaled"); hold on; title('337 Data'); 
    colorbar
           
end
figure;  image(Ratio,"CDataMapping","scaled"); hold on; title('Ratio Data'); colormap('jet'); colorbar;caxis([0,0.4]);
figure;  image(EoNData,"CDataMapping","scaled"); hold on; title('EoN Data'); 
figure;  image(A,"CDataMapping","scaled"); hold on;

% imcrop([352.5100    0.5100  671.9800  663.9800]);
colormap('jet'); a=colorbar; clim([0,1]);ylabel('Distance in y Direction / mm'); xlabel('Distance in y Direction / mm')
a.Label.String = 'Relative Intensity / a.u.';

% for nn = 1:Nfiles391-3
%     EoNPlotData = squeeze(EoNDataAll(nn,:,:));
%     figure;  image(EoNPlotData,"CDataMapping","scaled"); hold on; title('EoN Data'); 
%     colorbar, caxis([100,800])
% end
%% Create Gif of All the Data
CurrentTime = 0;
if LogicGif == 1
    for nn = 1:Nfiles337-2
         h = figure;
            %subplot(2,1,1)
            % A = squeeze(([OzoneAll(nn,:,:)+OzoneAll(nn-1,:,:)+OzoneAll(nn-2,:,:)] /3 ));
            EoNPlotData = squeeze(EoNDataAll(nn,:,:));
            %image(EoNPlotData,"CDataMapping","scaled"); hold on; title('EoN Data');
            figure;
            image(Ratio,"CDataMapping","scaled"); hold on; title('EoN Data');
            colorbar;
            colormap('jet')
            caxis([100,10000])
           
           
            colorbar
            ax = gca;
            ax.FontSize = 14;
            a = colorbar;
            ylabel(a,'Reduced Electric Field / Td','FontSize',14)
            hold on; % Don't let boundaries blow away the displayed image.



                CurrentTime = CurrentTime + 10;



                text(10,800,['Plasma On time = ',num2str(CurrentTime), ' \mu s'],'Color','white','FontSize',12)

            hold off;
            subplot(2,1,2)
            plot(VCData(:,1)*1000,VCData(:,3)*6)
            ylabel('Current / mA')
            xlabel('Time / ms')
            hold on 
            rectangle('Position',[CurrentTime/1000,-100,10/1000,200])
            axis([-30/1000,550/1000,-30,30])
            ax = gca;
            ax.FontSize = 14;
            ImageFilename = ['C:\Users\dptrosan\Documents\MATLAB\Results\ICCDKaptonEoN','Image_' , num2str(nn),'.png'];
            exportgraphics(h,ImageFilename)
            

    end
end



% %% Plot Last data Set
% if LogicReguclarPlotting ==1
%     figure
%     image(x,y,A,"CDataMapping","scaled"); hold on
%     xline(x(HighIndex))
%     xline(x(LowIndex))
%     colorbar
%     title('Absorbance')
%     xlabel('Distance in x Direction / mm')
%     ylabel('Distance in y Direction / mm')
%     ax = gca;
%     ax.FontSize = 24;
%     a = colorbar;
%     ylabel(a,'Absorbance / a.u.','FontSize',24)
% 
%     figure
%     image(x,y,Transmittance*100,"CDataMapping","scaled");
%     colorbar
%     title('Transmittance')
%     xlabel('Distance in x Direction / pixels')
%     ylabel('Distance in y Direction / pixels')
%     ax = gca;
%     ax.FontSize = 24;
%     a = colorbar;
%     ylabel(a,'Transmittance %','FontSize',24)
%     caxis([0,100])
% 
% 
%     figure
%     image(x,y,LightImage,"CDataMapping","scaled");
%     colorbar
%     title('Raw Data')
%     xlabel('Distance in x Direction / mm')
%     ylabel('Distance in y Direction / mm')
%     ax = gca;
%     ax.FontSize = 24;
%     a = colorbar;
%     ylabel(a,'Intensity / a.u.','FontSize',24)
%     hold on
%     for k = 1 : numberOfBoundaries
%         thisBoundary = boundaries{k}; % Get boundary for this specific blob.
%         x1 = thisBoundary(:,2)*ConversionFactor; % Column 2 is the columns, which is x.
%         y1 = thisBoundary(:,1)*ConversionFactor; % Column 1 is the rows, which is y.
%         plot(x1, y1, 'y-', 'LineWidth', 1); % Plot boundary in red.
%     end
% 
% 
%     figure
%     image(x,y,OzoneDensity,"CDataMapping","scaled"); hold on
%     xline(x(HighIndex))
%     xline(x(LowIndex))
%     colorbar
%     title('Ozone Density')
%     xlabel('Distance in x Direction / mm')
%     ylabel('Distance in y Direction / mm')
%     ax = gca;
%     ax.FontSize = 24;
%     a = colorbar;
%     ylabel(a,'Ozone / cm^{-3}','FontSize',24)
% 
%     %% AbleTransformation
%     % MeasurmentLine =  floor((AbelTransformLine)/ConversionFactor);
%     % DataSum = 0;
%     % for nn = 1:length()
%     % for i = 1:10
%     %     DataSum = OzoneDensity(MeasurmentLine+5 - i,:);
%     % end
%     % AbelData = DataSum/10;
% 
%     %% Plot Ozone All
%     i = 1;
%     TreatmentTime = [5,10,20,30,60];
% 
%     for nn = 1:9
%         A = squeeze(OzoneAll(nn,:,:));
%         if rem(nn,3) == 0
%             % Crop IMage to positions of intrest
%             cropnear = imcrop(A,rectnear);
%             cropfar = imcrop(A,rectfar);
% 
%             OzoneNear(i) = mean(cropnear,'All');
%             OzoneFar(i) = mean(cropfar,'all');
% 
%             figure
% 
%             image(x,y,A,"CDataMapping","scaled"); hold on
%             %xline(x(HighIndex))
%             %xline(x(LowIndex))
%             colorbar
%             title(['Ozone Density at 30 s and',num2str(TimeIndex(i)), ' Hertz'])
%             xlabel('Distance in x Direction / mm')
%             ylabel('Distance in y Direction / mm')
%             ax = gca;
%             ax.FontSize = 24;
%             a = colorbar;
%             caxis([1*10^12,2*10^14])
%             ylabel(a,'Ozone / cm^{-3}','FontSize',24)
%             hold on; % Don't let boundaries blow away the displayed image.
%             for k = 1 : numberOfBoundaries
%                 thisBoundary = boundaries{k}; % Get boundary for this specific blob.
%                 x1 = thisBoundary(:,2)*ConversionFactor; % Column 2 is the columns, which is x.
%                 y1 = thisBoundary(:,1)*ConversionFactor; % Column 1 is the rows, which is y.
%                 plot(x1, y1, 'y-', 'LineWidth', 1); % Plot boundary in red.
%             end
%             hold off;
%             OzoneInt(i) = sum(A,'all');
%             disp(['At t = ',num2str(TreatmentTime(i)), 'The integrated Ozone Number is ',num2str(OzoneInt(i)), 'Particles'])
%             i = i+1;
%         end
%     end
% 
%     figure;
%     plot(TimeIndex,OzoneNear); hold on
%     plot(TimeIndex,OzoneFar);
%     legend('Ozone Near Electrode','Ozone Far from Electrode')
%     ylabel('Ozone density / cm^{-3}')
%     xlabel('Frequency / Hz')
%     cat
%     % Sub Plot of Ozone
%     figure(1)
%     subplot(2,1,1)
%     image(x,y,squeeze(OzoneAll(15,:,:)),"CDataMapping","scaled"); hold on
%     %xline(x(HighIndex))
%     %xline(x(LowIndex))
%     colorbar
%     title('Ozone Density at 10 Seconds')
%     xlabel('Distance in x Direction / mm')
%     ylabel('Distance in y Direction / mm')
%     ax = gca;
%     ax.FontSize = 24;
%     a = colorbar;
%     caxis([10^14,5*10^16])
%     ylabel(a,'Ozone / cm^{-3}','FontSize',24)
% 
% 
%     subplot(2,1,2)
%     image(x,y,squeeze(OzoneAll(24,:,:)),"CDataMapping","scaled"); hold on
%     %xline(x(HighIndex))
%     %xline(x(LowIndex))
%     colorbar
%     title('Ozone Density at 60 Seconds')
%     xlabel('Distance in x Direction / mm')
%     ylabel('Distance in y Direction / mm')
%     ax = gca;
%     ax.FontSize = 24;
%     a = colorbar;
%     caxis([10^14,5*10^16])
%     ylabel(a,'Ozone / cm^{-3}','FontSize',24)
% 
% end
% %% Gif Plotting
% i = 1;
% CurrentTime = 0;




function Data = GetData
Data = [50	2.17631E-17	3.94892E-18	1.56615E-21	5.71144E-22	7.12965E-5	2.02272E-5
69.39	1.27334E-16	1.72921E-17	5.52049E-20	1.51681E-20	4.33291E-4	1.03591E-4
88.78	3.25572E-16	3.79094E-17	4.18525E-19	9.61254E-20	0.00129	2.7883E-4
108.2	5.78794E-16	6.19176E-17	1.57258E-18	3.1851E-19	0.00273	5.53315E-4
127.6	8.51512E-16	8.73944E-17	4.04082E-18	7.53782E-19	0.00477	9.24735E-4
146.9	1.12533E-15	1.12407E-16	8.27431E-18	1.44259E-18	0.0074	0.00139
166.3	1.38838E-15	1.3623E-16	1.45845E-17	2.39075E-18	0.01058	0.00193
185.7	1.6322E-15	1.58824E-16	2.31228E-17	3.62894E-18	0.01428	0.00256
205.1	1.85577E-15	1.79756E-16	3.39609E-17	5.14466E-18	0.01845	0.00325
224.5	2.06135E-15	1.99643E-16	4.71364E-17	6.93393E-18	0.02306	0.00401
243.9	2.24779E-15	2.17947E-16	6.25632E-17	8.97451E-18	0.02808	0.00482
263.3	2.41621E-15	2.34586E-16	8.01379E-17	1.12467E-17	0.03346	0.00569
282.7	2.56787E-15	2.49795E-16	9.97312E-17	1.37253E-17	0.03919	0.00661
302	2.7045E-15	2.63731E-16	1.21223E-16	1.6417E-17	0.04524	0.00758
321.4	2.82778E-15	2.76478E-16	1.44482E-16	1.92887E-17	0.05157	0.00858
340.8	2.93887E-15	2.88041E-16	1.69375E-16	2.23242E-17	0.05818	0.00962
360.2	3.03886E-15	2.98577E-16	1.9576E-16	2.55132E-17	0.06503	0.0107
379.6	3.12898E-15	3.0817E-16	2.23526E-16	2.88486E-17	0.07212	0.01182
399	3.21024E-15	3.16963E-16	2.52553E-16	3.23005E-17	0.07943	0.01296
418.4	3.28343E-15	3.24866E-16	2.82725E-16	3.58576E-17	0.08695	0.01414
437.8	3.34928E-15	3.32024E-16	3.13936E-16	3.95207E-17	0.09465	0.01534
457.1	3.40856E-15	3.3855E-16	3.46098E-16	4.32847E-17	0.10254	0.01657
476.5	3.462E-15	3.44481E-16	3.79139E-16	4.7133E-17	0.1106	0.01783
495.9	3.51002E-15	3.49835E-16	4.12961E-16	5.10422E-17	0.11882	0.01911
515.3	3.55306E-15	3.54682E-16	4.47479E-16	5.50131E-17	0.1272	0.02041
534.7	3.59151E-15	3.59009E-16	4.8263E-16	5.90493E-17	0.13573	0.02174
554.1	3.62596E-15	3.62992E-16	5.18365E-16	6.31467E-17	0.1444	0.02309
573.5	3.65669E-15	3.66568E-16	5.54625E-16	6.72888E-17	0.15321	0.02446
592.9	3.68398E-15	3.6975E-16	5.9135E-16	7.14535E-17	0.16215	0.02585
612.2	3.70817E-15	3.72648E-16	6.28499E-16	7.5663E-17	0.17122	0.02727
631.6	3.72947E-15	3.75212E-16	6.66014E-16	7.98984E-17	0.18041	0.0287
651	3.74812E-15	3.77535E-16	7.0387E-16	8.41688E-17	0.18972	0.03015
670.4	3.76443E-15	3.79613E-16	7.42034E-16	8.8487E-17	0.19915	0.03162
689.8	3.77857E-15	3.81448E-16	7.80475E-16	9.28116E-17	0.20869	0.03311
709.2	3.79076E-15	3.83054E-16	8.19164E-16	9.71647E-17	0.21834	0.03462
728.6	3.80109E-15	3.84517E-16	8.58077E-16	1.01525E-16	0.22809	0.03614
748	3.80974E-15	3.85784E-16	8.97174E-16	1.05903E-16	0.23795	0.03769
767.3	3.81673E-15	3.86831E-16	9.36398E-16	1.1025E-16	0.24791	0.03924
786.7	3.82216E-15	3.87731E-16	9.75776E-16	1.14629E-16	0.25798	0.04082
806.1	3.82613E-15	3.88479E-16	1.01526E-15	1.19004E-16	0.26814	0.04242
825.5	3.8288E-15	3.89074E-16	1.0548E-15	1.23403E-16	0.2784	0.04403
844.9	3.83027E-15	3.89544E-16	1.09444E-15	1.27811E-16	0.28876	0.04566
864.3	3.83069E-15	3.89924E-16	1.13415E-15	1.32215E-16	0.29922	0.0473
883.7	3.83009E-15	3.902E-16	1.17388E-15	1.36622E-16	0.30975	0.04897
903.1	3.82849E-15	3.90359E-16	1.21366E-15	1.41005E-16	0.32039	0.05064
922.4	3.82605E-15	3.90378E-16	1.25345E-15	1.45387E-16	0.33112	0.05233
941.8	3.8228E-15	3.90375E-16	1.29324E-15	1.49768E-16	0.34193	0.05404
961.2	3.81874E-15	3.90261E-16	1.33304E-15	1.54153E-16	0.35284	0.05577
980.6	3.814E-15	3.90093E-16	1.37281E-15	1.58484E-16	0.36383	0.0575
1000	3.80861E-15	3.89851E-16	1.41259E-15	1.62854E-16	0.37491	0.05925];
end

