%% PerimeterCoverage.m
% Computes, for each ICCD frame of the full electrode:
%   (1) Total plasma area [mm^2]
%   (2) Ignited open-perimeter fraction zeta(t) and length [mm]
% Then aligns zeta with the applied voltage so it can be compared with dQ/dV.
%
% Assumes the same Pimax CSV format / rotation / crop treatment as the
% single-hex script (ImageTreatment.mat). Runs at script level in the
% workspace like your other pipelines.

clear
close all

%% ---------------- User settings ----------------
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
BigHexDir = fullfile(RepoRoot,'Data','FullImagingDataSet','Fulldata','12_19_23Meanenergy','BigHex');
path            = fullfile(BigHexDir,'SweepOverTime','3372usSweep','Data');   % folder of frame CSVs
if ~exist(path,'dir'), path = uigetdir(pwd,'Choose the frame CSV folder'); end
BackgroundImagePath = fullfile(BigHexDir,'SweepOverTime','3372usSweep','BigHex2ussweep337LineBG 2023 December 19 17_28_14.csv');
ScopeFile = fullfile(RepoRoot,'Data','ICCD_ImagingData','ExampleTraces','BigHexOscope','PlasmaOn','bighex5v3000hz_ALL.csv');
% Gate delay -> scope time. The CSVs do not store the sweep's starting delay.
% 'auto' puts the end of the first emission hump on a voltage extremum of
% FirstHumpPolarity (discharge stops when |V| stops rising). The old fixed
% 1.71628 us offset leaves the light running ~50 us past each extremum.
ScopeTimeShift    = 'auto';      % 'auto' or a number in s (scope t = frame time + shift)
FirstHumpPolarity = 'negative';  % paper Fig 8
Fdrive = 3000;                   % Hz
AmpGain = 500;                   % V/V
PictureDt       = 2;         % us per frame
IgnitionCountsNegative  = 3000;       % counts above background = "plasma on"
IgnitionCountsPositive  = 1580;                             % (set from dark-frame noise: ~5-6 sigma)
BandWidthPx     = 4;         % px tolerance between emission ring and skeleton
MinBlobPx       = 10;        % despeckle: ignore blobs smaller than this
RefFrameMode    = 'LightImage';   % 'union' of all frames, or an integer frame no.
                             % Use a frame like your panel clsE (all cells lit,
                             % none saturated) if you want a specific one.
SkeletonMaskMode = 'refined'; % 'refined' (hexagons refitted to the light by
                              % ICCDFilterRatioEoN.m), 'load' (saved stamps, which sit
                              % ~3.6 px left of the real hexagons), 'stamp', 'manual',
                              % 'interpolate'
PlotAll = 'all';             % Plot all the data
LightImagePath = fullfile(BigHexDir,'LightImageBigHex 2023 December 19 17_31_04.csv');
ImageAdjustments = 'no';


%% ---------------- Read background ----------------
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [2, Inf];  opts.Delimiter = ",";
opts.VariableNames = ["Y","X","I"];
opts.VariableTypes = ["double","double","double"];
opts.ExtraColumnsRule = "ignore"; opts.EmptyLineRule = "read";

BackgroundMatrix = table2array(readtable(BackgroundImagePath, opts));
Background = reshape(BackgroundMatrix(:,3), 1024, 1024)';   % row-major restructure
% Background = imrotate(Background, RotationDegrees, 'bilinear', 'crop');

T = table2array(readtable(LightImagePath, opts));
LightImage = reshape(T(:,3), 1024, 1024)';

%% ---------------- Create Conversion Factor ---------------------------
%% Find ConversionFactor if Logic Conversion Factor is set to 'yes'
if strcmp(ImageAdjustments,'yes') == 1
    % Rotate Image
    figure;
    imshow(LightImage,[1,max(max(LightImage))])
    LightImageTemp = LightImage;
    % Crop Image

    disp('Do you want to skip Image Cropping')
    LogicYes = input('Yes or no','s');
    if strcmpi(LogicYes,'no') == 1
        clf
        figure;
        image(LightImageTemp,"CDataMapping","scaled");colormap(gray)
        h = imrect(gca,[500 500 500 500]);
        Position = wait(h); % returns coordinates in "position" when user doubleclicks on rectangle
        CropLightImage = imcrop(LightImageTemp,Position);
        % [CropLightImage,croprect] = imcrop;
    else
        CropLightImage = LightImageTemp;
        Position = [1,1,1023,1023];
    end

    %% Conversion Factor
    clf
    figure;
    image(CropLightImage,"CDataMapping","scaled");colormap(gray)
    % Final Image Adjustments
    for n = 1:5
        disp(['Please indentify two distance points ',num2str(5-(n-1)),' more times'])
        [x,y] = ginput(2);
        Pixels(n) = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
        disp(['Pixel length = ',num2str(Pixels)])
    end
    Distance = input('Distance in mm = ');
    ConversionFactor = mean(Distance./Pixels); %% Conversion factor in mm/pixel
    ImageHandeling = struct;
    ImageHandeling.Crop = Position;
    ImageHandeling.ConversionFactor = ConversionFactor;
    save(fullfile(RepoRoot,'BiggHexImaging','ImageTreatment.mat'),"ImageHandeling")
    CropRect = Position;
else
    load(fullfile(RepoRoot,'BiggHexImaging','ImageTreatment.mat'))   % same rotate/crop/ConversionFactor as before
    CropRect         = ImageHandeling.Crop;
    ConversionFactor = ImageHandeling.ConversionFactor;  % mm/px

end
Background = imcrop(Background, CropRect);


%% ---------------- Read all frames ----------------
FileList = dir(fullfile(path,'*.csv'));
Nfiles   = length(FileList);
[Ny,Nx]  = size(Background);
FinalData = zeros(Nfiles, Ny, Nx);

for nn = 1:Nfiles
    T = table2array(readtable(fullfile(FileList(nn).folder,FileList(nn).name), opts));
    Img = reshape(T(:,3), 1024, 1024)';
    % Img = imrotate(Img, RotationDegrees, 'bilinear', 'crop');
    Img = imcrop(Img, CropRect);
    FinalData(nn,:,:) = medfilt2(Img - Background, [3,3]);
end
TimeIndex = (0:Nfiles-1) * PictureDt;   % us

%% ---------------- Build the open-perimeter skeleton ----------------
% Reference image: every cell edge lit at some point in the cycle.
if strcmp(RefFrameMode,'union') == 1
    RefImage = squeeze(max(FinalData,[],1));
elseif strcmp(RefFrameMode,'LightImage') == 1
    LightImage = imcrop(LightImage,CropRect);
    RefImage = LightImage;
else
    RefImage = squeeze(FinalData(RefFrameMode,:,:));
end

if strcmp(SkeletonMaskMode,'manual') == 1
    figure; imshow(RefImage, [0, max(RefImage(:))]); hold on
    title('Trace each hexagon perimeter. Double-click / Enter to finish each one.')

    PerimMask = false(size(RefImage));   % accumulator

    % ---- Option A: polygon per hexagon (click the 6 corners) ----
    Ncells = 33;                         % however many you have
    for k = 1:Ncells
        h = drawpolygon('Color','c');    % click 6 vertices, double-click to close
        ringMask = bwperim(createMask(h));   % OUTLINE only, not filled
        PerimMask = PerimMask | ringMask;
    end


    % ---- Option B: freehand trace instead (swap into the loop) ----
    % h = drawfreehand('Color','c','Closed',true);
    % ringMask = bwperim(createMask(h));
    % PerimMask = PerimMask | ringMask;

    Skeleton = PerimMask;
    PerimTotalPx = nnz(Skeleton);
    save(fullfile(RepoRoot,'BiggHexImaging','PerimeterSkeleton.mat'),'Skeleton','PerimTotalPx')

elseif strcmp(SkeletonMaskMode,'interpolate')
    RefMask = RefImage > IgnitionCountsNegative;
    RefMask = bwareaopen(RefMask, MinBlobPx);
    RefMask = imclose(RefMask, strel('disk',2));
    Skeleton = bwskel(RefMask);              % 1-px-wide trace of every hex edge


    % Remove the outer carrier arc / anything that is not a hexagon ring:
    % keep connected skeleton pieces whose bounding box matches a single cell.
    CC = bwconncomp(Skeleton);
    S  = regionprops(CC,'BoundingBox');
    CellSizePx = 2*4.56 / ConversionFactor;  % hex cell width, side 4.56 mm (paper)
    Keep = false(1,CC.NumObjects);
    for k = 1:CC.NumObjects
        bb = S(k).BoundingBox;
        Keep(k) = bb(3) < 2*CellSizePx && bb(4) < 2*CellSizePx;
    end
    Skeleton = ismember(labelmatrix(CC), find(Keep));

    % ---- Stamp hexagons and manually select centers
    %% Manual perimeter by center-tapping — stamps a regular hexagon per click
elseif strcmp(SkeletonMaskMode,'stamp') == 1
    disp('Stamp Mode Activated')
    SideMM   = 4;
    R_px     = SideMM / ConversionFactor;
    Orient   = 'flat';
    if strcmp(Orient,'flat')
        ang = (0:5)*60;
    else
        ang = (0:5)*60 + 30;
    end
    ang = deg2rad(ang);
    fig = figure;
    imshow(RefImage, [0, max(RefImage(:))]); hold on
    title("Click each cell CENTER. Press 'q' to stop.")

    PerimMask = false(size(RefImage));
    Rings   = {};
    Centers = [];
    k = 0;

    % ---------- Phase 1: stamp regular hexagons ----------
    while true
        [xc, yc, button] = ginput(1);
        if isempty(button) || button == double('q')
            break
        end
        vx = xc + R_px*cos(ang);
        vy = yc + R_px*sin(ang);
        hexMask  = poly2mask(vx, vy, size(RefImage,1), size(RefImage,2));
        ringMask = bwperim(hexMask);
        k = k + 1;
        Rings{k}     = ringMask;
        Centers(k,:) = [xc, yc];
        PerimMask    = PerimMask | ringMask;
        plot([vx vx(1)],[vy vy(1)],'c-','LineWidth',1);
        title(sprintf('%d cells placed. Press q to stop stamping.', k));
    end
    fprintf('Stamped %d hexagons.\n', k);

    % ---------- Phase 2: freehand additions ----------
    title("Freehand: trace extra perimeter. Enter after each. Press 'q' to finish.")
    FreeRings = {};
    m = 0;
    while true
        h = drawfreehand('Color','y','Closed',false);   % trace a segment/loop
        if isempty(h) || isempty(h.Position) || size(h.Position,1) < 2
            break                                        % empty draw = done
        end
        % outline of the freehand region (use bwperim for closed shapes;
        % for open traces the mask itself is the thin band)
        fMask = createMask(h);
        if h.Closed
            fMask = bwperim(fMask);
        end
        m = m + 1;
        FreeRings{m} = fMask;
        PerimMask = PerimMask | fMask;
        title(sprintf('%d freehand traces added. Press q to finish.', m));
        % ask whether to continue
        btn = questdlg('Add another freehand trace?','Continue','Yes','No','No');
        if strcmp(btn,'No'), break, end
    end
    fprintf('Added %d freehand traces.\n', m);

    % ---------- Finalize (after BOTH phases) ----------
    Skeleton     = PerimMask;
    PerimTotalPx = nnz(Skeleton);
    fprintf('Total: %d hexagons + %d freehand, %d perimeter px\n', k, m, PerimTotalPx);
    save(fullfile(RepoRoot,'BiggHexImaging','PerimeterSkeleton.mat'), ...
         'Skeleton','PerimTotalPx','Rings','Centers','FreeRings','SideMM')
% --------- Load previous saved skeleton
elseif strcmp(SkeletonMaskMode,'load')
    load(fullfile(RepoRoot,'BiggHexImaging','PerimeterSkeleton.mat'))
elseif strcmp(SkeletonMaskMode,'refined')
    RefFile = fullfile(RepoRoot,'Data','FullImagingDataSet','Results','EoN_BigHex','HexagonGeometryRefined.mat');
    if ~exist(RefFile,'file')
        error('%s not found: run ImagingAnalysis/ICCDFilterRatioEoN.m (LogicSides = 1) first.',RefFile);
    end
    R = load(RefFile,'SkelAll');
    Skeleton = R.SkelAll;
end




PerimTotalPx = nnz(Skeleton);
PerimTotalMM = PerimTotalPx * ConversionFactor;
disp(['Total open perimeter = ', num2str(PerimTotalMM), ' mm'])

% Sanity check against geometry: N_cells * 6 * side length should be close.

%% ---------------- Per-frame coverage and area ----------------
Band = strel('disk', BandWidthPx);
Zeta        = zeros(1,Nfiles);   % ignited perimeter fraction
PerimLitMM  = zeros(1,Nfiles);   % ignited perimeter length, mm
PlasmaArea  = zeros(1,Nfiles);   % total plasma area, mm^2

for nn = 1:Nfiles
    Img  = squeeze(FinalData(nn,:,:));
    if TimeIndex(nn) < 166
        IgnitionCounts = IgnitionCountsNegative;
    else
        IgnitionCounts = IgnitionCountsPositive;
    end

    Mask = bwareaopen(Img > IgnitionCounts, MinBlobPx);
    PlasmaArea(nn) = nnz(Mask) * ConversionFactor^2;
    TotalInt(nn) = sum(sum(Img(Mask)));
    LitSkel = imdilate(Mask, Band) & Skeleton;
    Zeta(nn)       = nnz(LitSkel) / PerimTotalPx;
    PerimLitMM(nn) = nnz(LitSkel) * ConversionFactor;

end

%% ---------------- QC overlay (spot-check a few frames) ----------------
CheckFrames = round(linspace(2, Nfiles-1, 6));
figure
for k = 1:6
    nn = CheckFrames(k);
    subplot(2,3,k)
    Img  = squeeze(FinalData(nn,:,:));
    imshow(Img, [0, 4000]); colormap(gca,'gray'); hold on

    Mask = bwareaopen(Img > IgnitionCounts, MinBlobPx);
    Lit  = imdilate(Mask, Band) & Skeleton;
    Dark = Skeleton & ~Lit;
    [yy,xx] = find(Lit);  plot(xx,yy,'.g','MarkerSize',2)
    [yy,xx] = find(Dark); plot(xx,yy,'.r','MarkerSize',2)
    title(['t = ', num2str(TimeIndex(nn)), ' \mus,  \zeta = ', num2str(Zeta(nn),'%.2f')])
end

% Plot all the figures
if strcmp(PlotAll,'all') == 1
for nn = 1:Nfiles
    figure
    Img  = squeeze(FinalData(nn,:,:));
    imshow(Img, [0, 4000]); colormap(gca,'gray'); hold on
    Mask = bwareaopen(Img > IgnitionCounts, MinBlobPx);
    Lit  = imdilate(Mask, Band) & Skeleton;
    Dark = Skeleton & ~Lit;
    [yy,xx] = find(Lit);  plot(xx,yy,'.g','MarkerSize',2)
    [yy,xx] = find(Dark); plot(xx,yy,'.r','MarkerSize',2)
    title(['t = ', num2str(TimeIndex(nn)), ' \mus,  \zeta = ', num2str(Zeta(nn),'%.2f')])
    
end

% Mask Images 
for nn = 1:Nfiles
    figure
    Img  = squeeze(FinalData(nn,:,:));
    Mask = bwareaopen(Img > IgnitionCounts, MinBlobPx);
    imshow(Mask, [0, 1]); colormap(gca,'gray'); hold on
    Lit  = imdilate(Mask, Band) & Skeleton;
    Dark = Skeleton & ~Lit;
    title(['t = ', num2str(TimeIndex(nn)), ' \mus,  \zeta = ', num2str(Zeta(nn),'%.2f')])
    
end


% Light Image 
figure
imshow(RefImage, [0, max(RefImage(:))]); hold on
[yy, xx] = find(Skeleton);
plot(xx, yy, '.', 'Color', [0 1 0], 'MarkerSize', 3);   % green skeleton dots
title('Light image with traced perimeter')

end

%% ---------------- Zeta vs time / voltage ----------------
figure
yyaxis left
plot(TimeIndex, Zeta, 'LineWidth', 2)
ylabel('Ignited perimeter fraction \zeta')
yyaxis right
plot(TimeIndex, PlasmaArea, 'LineWidth', 2)
ylabel('Plasma area / mm^2')
xlabel('Time / \mus')
ax = gca; ax.FontSize = 18;

% Align with the oscilloscope trace
Scope = readmatrix(ScopeFile,'NumHeaderLines',11);
t = Scope(:,1);
v = Scope(:,2)*AmpGain;
dt = median(diff(t));
vs = movmean(v,round(5e-6/dt));
[~,iMax] = findpeaks( vs,'MinPeakDistance',round(0.6/Fdrive/dt));
[~,iMin] = findpeaks(-vs,'MinPeakDistance',round(0.6/Fdrive/dt));
if ischar(ScopeTimeShift) && strcmpi(ScopeTimeShift,'auto')
    RawTot = squeeze(sum(FinalData,[2 3]))';   % all light, not the thresholded mask
    On = RawTot > 0.01*max(RawTot);
    HumpEnd = find(On & ~[On(2:end) false],1,'first');
    if strcmpi(FirstHumpPolarity,'negative'), tExt = t(iMin); else, tExt = t(iMax); end
    Cand = tExt - (TimeIndex(HumpEnd) + PictureDt)*1e-6;
    Cand = Cand(Cand + TimeIndex(1)*1e-6 >= t(1) & Cand + (TimeIndex(end)+PictureDt)*1e-6 <= t(end));
    [~,j] = min(abs(Cand));
    ScopeTimeShift = Cand(j);
end
fprintf('Scope time shift = %.1f us\n',ScopeTimeShift*1e6);
FrameTimes = (TimeIndex + PictureDt/2)*1e-6 + ScopeTimeShift;   % s, gate centre
Vframe = interp1(t, v, FrameTimes, 'linear', NaN);
VframeOld = interp1(t, v, TimeIndex*1e-6 - 1.71628e-6, 'linear', NaN);   % old offset, for comparison

figure
plot(Vframe, Zeta, 'o-', 'LineWidth', 2)
xlabel('Applied Voltage / V'); ylabel('\zeta')
ax = gca; ax.FontSize = 18;
% If ignition thresholds along the perimeter are ~Gaussian distributed,
% zeta(V) on the rising quarter-cycle should follow a cumulative normal:
%   zeta(V) = 0.5*(1 + erf((V - mu)/(sqrt(2)*sigmaV)))
% Fit mu, sigmaV and compare with the width of the dQ/dV transition.
figure; 
plot(TimeIndex, TotalInt)


OriginZeta = transpose([TimeIndex; Vframe; Zeta; PerimLitMM; PlasmaArea]);


% Scalled Total Intensity 
for n = 1:Nfiles
    if TimeIndex(n) < 167 
        ScaledFactor = 1.96; 
        A(n) = TotalInt(n)/ScaledFactor;
    else 
                ScaledFactor = 1; 
        A(n) = TotalInt(n)/ScaledFactor;
    end
end
figure;
plot(TimeIndex,A)
   