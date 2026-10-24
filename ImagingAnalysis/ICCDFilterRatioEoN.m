%% ICCD Filter Ratio Imaging (337 / 391 nm)
%
% Repo-relative rework of ICCdFilterAnalsysisCode.m for the Big Hex 2 us
% gate sweeps of 12/19/23 (3000 Hz, 2.5 kV, 3000 accumulations). For every
% gate position it builds the N2(C)/N2+(B) emission ratio map (337 nm band
% photons over 391 nm N2+ band photons, filter corrected and anchored to the
% echelle spectrum), and optionally:
%   ReportEoN    - also converts the ratio to E/N with the calibration table
%   LogicGif     - writes the ratio map + current trace GIF (one frame per gate)
%   LogicDelayed - compares the 391/337 ratio of hexagons that ignite late
%                  with the hexagons that ignite first, per half-wave
%
% CHANGES FROM ICCdFilterAnalsysisCode.m
%   - Filter factors come from the echelle spectrum (EchelleFilterFactors.m)
%     instead of the peak transmissions. FilterFactorMode = 'legacy' restores
%     0.2815 / 0.4134 / 0.4533.
%   - SPS leakage into the 391 filter is subtracted as Kappa*S337 instead of
%     multiplying S391 by a fixed N2+ fraction (ContaminationMode). The fixed
%     fraction is only right at the E/N of the spectrum it came from.
%   - The 391 stack is registered to the 337 stack (the filter swap shifts the
%     image by a few pixels; a pixel-wise ratio is very sensitive to that).
%   - Same crop / conversion factor as the perimeter code (ImageTreatment.mat)
%     so the hexagon skeleton (PerimeterSkeleton.mat) lines up.
%   - The GIF loop plotted the Ratio image under an E/N colorbar, opened two
%     figures per frame and exported the empty one. Fixed.
%   - Echelle anchoring (AnchorToEchelle): the 391 signal is scaled once so the
%     electrode-integrated, full-period raw 391/337 equals the echelle's. Both
%     instruments see the same time-integrated light, so this one factor
%     absorbs accumulation, gain and QE differences between the two sweeps and
%     the echelle, and the echelle's N2+ share and Kappa then apply directly.
%   - The reported quantity is the N2(C)/N2+(B) emission ratio, which needs no
%     kinetic model. It falls as the field rises.
%   - E/N calibration (ReportEoN = 1; GetData, 4th order fit of column 1 vs
%     column 6) is unchanged. NOTE: column 6 is ~k(N2+ B)/k(N2 C), so the mapping treats
%     the corrected photon ratio as the excitation-rate ratio, with no
%     lifetime/quenching/branching correction. Absolute E/N carries that
%     assumption; relative comparisons (e.g. delayed vs early hexagons) do not
%     depend on it because the mapping is monotonic.

clear
% No 'close all': figures from this run stay open (and every figure is also
% saved as .fig next to its .png, so it can be reopened and edited later).

%% Section 1 - Options
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
DataRoot = fullfile(RepoRoot,'Data','FullImagingDataSet','Fulldata','12_19_23Meanenergy','BigHex');
SweepDir337 = fullfile(DataRoot,'SweepOverTime','3372usSweep');
SweepDir391 = fullfile(DataRoot,'SweepOverTime','3912usSweep');
ScopeFile   = fullfile(RepoRoot,'Data','ICCD_ImagingData','ExampleTraces','BigHexOscope','PlasmaOn','bighex5v3000hz_ALL.csv');
SpecFile    = fullfile(RepoRoot,'Data','ICCD_ImagingData','Echelledata','AbsoultelyCalibratedtrial', ...
                       'Absoultecalibrationtest','Test','KaptonSmallHex.spa');
TransFile   = fullfile(RepoRoot,'Data','Stuff','FilterTransmitance.csv');
TreatFile   = fullfile(RepoRoot,'BiggHexImaging','ImageTreatment.mat');
SkelFile    = fullfile(RepoRoot,'BiggHexImaging','PerimeterSkeleton.mat');
% All results go to Output/, one folder per kind of result
OutRoot = fullfile(RepoRoot,'Output');
OutDir.Cache      = fullfile(OutRoot,'Cache');         % stacks, derived data, geometry
OutDir.Ratio      = fullfile(OutRoot,'EmissionRatio'); % GIF, per-gate emission and ratio
OutDir.Sides      = fullfile(OutRoot,'SideAnalysis');  % per-side rectangles (Section 10)
OutDir.Perimeter  = fullfile(OutRoot,'Perimeter');     % open perimeter ratio
OutDir.Superseded = fullfile(OutRoot,'Superseded');    % Sections 8 and 9
CacheFile   = fullfile(OutDir.Cache,'BigHexFilterStacks.mat');

UseCache          = 1;          % reuse the cropped, background-subtracted stacks
FilterFactorMode  = 'echelle';  % 'echelle' or 'legacy'
ContaminationMode = 'fraction'; % 'fraction' (S391*FNSFraction, the established
                                % method) or 'subtract' (S391 - Kappa*S337).
                                % 'subtract' is the better leakage model; with
                                % AnchorToEchelle = 1 the echelle Kappa is in
                                % the right units, but it nearly cancels the
                                % negative half-wave's N2+ signal.
AnchorToEchelle   = 1;          % scale ICCD 391 so the full-period raw 391/337
                                % equals the echelle's (see header)
ReportEoN         = 0;          % 1 = also report E/N (Td) next to the ratio
PhotonUnits       = 1;          % echelle .spa in energy units -> photons
KappaOverride     = [];         % set a number to override the echelle Kappa

Fdrive        = 3000;           % Hz
PictureDt     = 2;              % us between gate positions
GateWidth     = 2;              % us
% Gate delay -> scope time. The CSVs do not record the sweep's starting delay,
% and the 1.7 us offset in PerimeterCodev2.m leaves the light running ~50 us
% past each voltage extremum, which a DBD cannot do. 'auto' lines the end of
% the first emission hump up with a voltage extremum of the polarity below.
ScopeTimeShift    = 'auto';     % 'auto' or a number in s (scope t = gate + shift)
FirstHumpPolarity = 'negative'; % paper Fig 8: first half-wave imaged is negative
AmpGain       = 500;            % V/V, Trek voltage monitor
CurrentGain   = 6;              % mA/V, Trek current monitor

RatioSmoothPx = 5;              % box filter (px) on both images before the ratio
Min337Ratio   = 1500;           % counts (after smoothing) to compute a ratio
Min391Ratio   = 300;            % counts (after smoothing) to compute a ratio
RatioMax      = 0.6;            % cap on N2+(B)/N2(C) (noise), same as the original code
CBLimits      = [2 30];         % N2(C)/N2+(B) colour scale (log)
EoNLimits     = [100 1000];     % Td, colour scale (ReportEoN = 1)

LogicGif      = 1;
GifFile       = fullfile(OutDir.Ratio,'BigHexEmissionRatio.gif');
GifQuantity   = 'ratio';        % 'ratio', or 'eon' (needs ReportEoN = 1)
GifDelay      = 0.12;           % s per GIF frame
GifFrames     = [];             % [] = every gate position

LogicDelayed  = 0;              % per-hexagon version; superseded by LogicSides and still
                                % uses the unrefined stamps (offset ~3.6 px)
IgnitionCountsNegative = 3000;  % 337 counts, same thresholds as PerimeterCodev2.m
IgnitionCountsPositive = 1580;
BandWidthPx   = 4;              % px tolerance around each traced hexagon edge
CellDilatePx  = 3;              % px added outside each hexagon for its cell mask
LitFraction   = 0.5;            % cell has ignited when this share of its edge is lit
DelayedFrames = 2;              % delayed = ignites >= this many gates after the MEDIAN cell
                                % (the Spearman correlation needs no threshold)
OnsetFrames   = 2;              % gates after each cell's own ignition for the onset ratio

LogicDelayedPixel   = 0;        % pixel-by-pixel version of the delayed-ignition test
                                % (off: plasma expansion makes outer pixels look late)
PixelOnCounts       = [3000 1580]; % [negative positive] smoothed 337 counts a pixel must
                                % reach at some point in the half-wave to count as "turns on"
PixelIgnFraction    = 0.25;     % pixel ignites at the first gate above this share of its
                                % own half-wave maximum (relative, so dim pixels are not
                                % called late just for being dim)
PixelDelayThreshold = 6;        % us; delayed = ignites this much after the reference
PixelRefPercentile  = 10;       % reference ignition = this percentile of all ignition times
RadiusBins          = 6;        % radial strata for the position-controlled comparison

LogicSides       = 1;           % per-side rectangles + open perimeter ratio
RefineHexagons   = 1;           % refit each stamped hexagon's centre and size to the light.
                                % The stamps in PerimeterSkeleton.mat sit ~3.6 px left of the
                                % real hexagons, which made the right-hand sides look late.
MaxRefinePx      = 8;           % reject a refit that moves an edge further than this
ExtentCounts     = [1000 1000]; % [negative positive] mean 337 counts per pixel row that
                                % still counts as plasma for the extent (Fig 3 threshold)
DerivedFile      = fullfile(OutDir.Cache,'BigHexDerived.mat');   % everything Fig13Compare.m needs
SideTrim         = 0.15;        % share of the side length dropped at each corner
SideDepthMM      = 2.0;         % rectangle depth into the hexagon opening, mm (extent ~1.6 mm)
SideOutMM        = 0.4;         % rectangle depth outward over the foil edge, mm
SideEdgeHalfPx   = 1;           % half-width of the edge line used for the lit fraction, px
SideLitFraction  = 0.5;         % side has ignited when this share of its edge line is lit
SideDelayThreshold = 4;         % us; delayed = ignites this much after the median side
SideMatchedPct   = 95;          % matched window starts once this % of sides have ignited

for f = fieldnames(OutDir)'
    if ~exist(OutDir.(f{1}),'dir'), mkdir(OutDir.(f{1})); end
end

%% Section 2 - Filter factors
switch lower(FilterFactorMode)
case 'echelle'
    FF = EchelleFilterFactors(SpecFile,TransFile,0,PhotonUnits);
    FilterFactor337     = FF.T337eff;
    FilterFactor391     = FF.T391eff;
    Kappa               = FF.Kappa;
    Factor391Percentage = FF.FNSFraction;
    EchelleRaw          = (FF.SPS391 + FF.FNS391)/FF.S337;   % full-period raw 391/337
case 'legacy'
    FilterFactor337     = 0.2815;
    FilterFactor391     = 0.4134;
    Factor391Percentage = 0.4533;
    Kappa               = NaN;
    EchelleRaw          = NaN;
otherwise
    error('FilterFactorMode must be ''echelle'' or ''legacy''.');
end
if ~isempty(KappaOverride), Kappa = KappaOverride; end
if strcmpi(ContaminationMode,'subtract') && isnan(Kappa)
    error('ContaminationMode ''subtract'' needs Kappa (use the echelle mode or KappaOverride).');
end

% E/N calibration, unchanged from the original code
ReducedElectricFieldData = GetData;
EoNCoeff = polyfit(ReducedElectricFieldData(:,6),ReducedElectricFieldData(:,1),4);
RatioRange = [min(ReducedElectricFieldData(:,6)) max(ReducedElectricFieldData(:,6))];

%% Section 3 - Load the image stacks (cropped, background subtracted, 3x3 median)
IT = load(TreatFile);
CropRect = IT.ImageHandeling.Crop;
ConversionFactor = IT.ImageHandeling.ConversionFactor;   % mm/px

if UseCache == 1 && exist(CacheFile,'file')
    Stk = load(CacheFile);
    fprintf('Loaded cached stacks from %s\n',CacheFile);
else
    Stk = struct;
    Lines = {'337',SweepDir337; '391',SweepDir391};
    for L = 1:2
        bgf = dir(fullfile(Lines{L,2},'*BG*.csv'));
        BG  = imcrop(ReadPimaxCSV(fullfile(bgf(1).folder,bgf(1).name)),CropRect);
        fl  = dir(fullfile(Lines{L,2},'Data','*.csv'));
        [~,ix] = sort({fl.name}); fl = fl(ix);
        S = zeros([size(BG) numel(fl)],'single');
        for n = 1:numel(fl)
            Img = imcrop(ReadPimaxCSV(fullfile(fl(n).folder,fl(n).name)),CropRect);
            S(:,:,n) = single(medfilt2(Img - BG,[3 3]));
        end
        Stk.(['S' Lines{L,1}])     = S;
        Stk.(['Files' Lines{L,1}]) = {fl.name};
        Stk.(['BG' Lines{L,1}])    = single(BG);
        fprintf('%s nm: %d frames loaded\n',Lines{L,1},numel(fl));
    end
    Stk.CropRect = CropRect;
    save(CacheFile,'-struct','Stk','-v7.3');
end
S337 = Stk.S337;     % single, to keep the full stacks in memory
S391 = Stk.S391;
[Ny,Nx,Nframes] = size(S337);
if size(S391,3) ~= Nframes
    error('337 and 391 sweeps have different frame counts (%d vs %d).',Nframes,size(S391,3));
end
GateDelay = (0:Nframes-1)*PictureDt;   % us

%% Section 4 - Register the 391 stack onto the 337 stack
% Translation only, from the time-summed images. The filters are swapped in
% front of the lens, which can shift the image by a few pixels.
Sum337 = sum(S337,3);
Sum391 = sum(S391,3);
tform = imregcorr(Sum391/max(Sum391(:)),Sum337/max(Sum337(:)),'translation');
Shift = tform.T(3,1:2);
fprintf('391 -> 337 registration shift: dx = %.2f px, dy = %.2f px\n',Shift(1),Shift(2));
Ref = imref2d([Ny Nx]);
for n = 1:Nframes
    S391(:,:,n) = imwarp(S391(:,:,n),tform,'OutputView',Ref,'FillValues',0);
end

%% Section 5 - Oscilloscope trace and gate timing
Scope = readmatrix(ScopeFile,'NumHeaderLines',11);
tScope   = Scope(:,1);                       % s
VScope   = Scope(:,2)*AmpGain;               % V
IScope   = Scope(:,3)*CurrentGain;           % mA
dtScope  = median(diff(tScope));
% Remove the current probe offset over whole periods only
NperT    = round(1/Fdrive/dtScope);
Nwhole   = floor(numel(IScope)/NperT)*NperT;
IScope   = IScope - mean(IScope(1:Nwhole));

% Emission humps: runs of gates with total 337 light above 1% of the maximum
Tot337 = double(squeeze(sum(S337,[1 2])));
On     = Tot337(:)' > 0.01*max(Tot337);
HumpEnd   = find(On & ~[On(2:end) false]);
HumpStart = find(On & ~[false On(1:end-1)]);
VSmooth = movmean(VScope,round(5e-6/dtScope));
[~,iMax] = findpeaks(VSmooth, 'MinPeakDistance',round(0.6/Fdrive/dtScope));
[~,iMin] = findpeaks(-VSmooth,'MinPeakDistance',round(0.6/Fdrive/dtScope));
if ischar(ScopeTimeShift) && strcmpi(ScopeTimeShift,'auto')
    % Discharge stops at the voltage extremum, so the end of hump 1 (+ one
    % gate width) sits on an extremum of FirstHumpPolarity. Of the matching
    % extrema, take the one that keeps the whole sweep inside the scope record.
    if strcmpi(FirstHumpPolarity,'negative'), tExt = tScope(iMin); else, tExt = tScope(iMax); end
    tEnd1 = (GateDelay(HumpEnd(1)) + GateWidth)*1e-6;
    Cand  = tExt - tEnd1;
    Fits  = Cand + GateDelay(1)*1e-6 >= tScope(1) & Cand + (GateDelay(end)+GateWidth)*1e-6 <= tScope(end);
    Cand  = Cand(Fits);
    [~,j] = min(abs(Cand));
    ScopeTimeShift = Cand(j);
end
GateStart = GateDelay*1e-6 + ScopeTimeShift;  % s, scope time of each gate opening
GateMid   = GateStart + GateWidth*1e-6/2;
VGate     = interp1(tScope,VScope,GateMid,'linear',NaN);
dVGate    = interp1(tScope,gradient(VSmooth,dtScope),GateMid,'linear',NaN);
% Negative half-wave = powered electrode voltage falling. Each emission hump is
% assigned as a whole, from the slope at its middle gate: its last gates sit
% on the voltage extremum where dV/dt ~ 0 and the sign there is noise.
HalfWave  = ones(1,Nframes);                 % +1 positive, -1 negative
HalfWave(dVGate < 0) = -1;
for h = 1:numel(HumpEnd)
    Mid = round((HumpStart(h) + HumpEnd(h))/2);
    HalfWave(HumpStart(h):HumpEnd(h)) = sign(dVGate(Mid));
end
fprintf('Scope time shift: %.1f us (scope t = gate delay + shift)\n',ScopeTimeShift*1e6);
for h = 1:numel(HumpEnd)
    tEnd = GateStart(HumpEnd(h)) + GateWidth*1e-6;
    [dMax,~] = min(abs(tScope(iMax) - tEnd));
    [dMin,~] = min(abs(tScope(iMin) - tEnd));
    if HalfWave(HumpEnd(h)) < 0, Pol = 'negative'; else, Pol = 'positive'; end
    fprintf('  Hump %d: gates %d-%d us, %s half-wave, ends %.1f us from the nearest V max, %.1f us from the nearest V min\n', ...
            h,GateDelay(HumpStart(h)),GateDelay(HumpEnd(h)),Pol,dMax*1e6,dMin*1e6);
end

%% Section 6 - Echelle anchor and ratio maps
Tot391 = double(squeeze(sum(S391,[1 2])));
if AnchorToEchelle == 1 && ~isnan(EchelleRaw)
    AnchorGain = EchelleRaw/(sum(Tot391(On))/sum(Tot337(On)));
else
    AnchorGain = 1;
end
fprintf('Echelle anchor: full-period raw 391/337 %.3f -> %.3f (391 x %.3f)\n', ...
        sum(Tot391(On))/sum(Tot337(On)),AnchorGain*sum(Tot391(On))/sum(Tot337(On)),AnchorGain);

if RatioSmoothPx > 1
    A = imboxfilt3(S337,[RatioSmoothPx RatioSmoothPx 1]);
    B = imboxfilt3(S391,[RatioSmoothPx RatioSmoothPx 1]);
else
    A = S337;  B = S391;
end
Valid391 = B >= Min391Ratio;   % threshold on the measured counts, before anchoring
B = AnchorGain*B;
switch lower(ContaminationMode)
case 'subtract'
    BFNS = B - Kappa*A;
    Bright = A > Min337Ratio;
    RawMin = prctile(B(Bright)./A(Bright),5);
    if Kappa > RawMin
        warning(['Kappa = %.3f is above the 5th percentile of the raw 391/337 ' ...
                 '(%.3f): the corrected N2+ signal goes negative there. Kappa ' ...
                 'needs to be in ICCD units.'],Kappa,RawMin);
    end
case 'fraction'
    BFNS = B*Factor391Percentage;
otherwise
    error('ContaminationMode must be ''subtract'' or ''fraction''.');
end
Ratio = (BFNS/FilterFactor391)./(A/FilterFactor337);     % N2+(B)/N2(C) photons
Valid = A >= Min337Ratio & Valid391 & Ratio > 0 & Ratio < RatioMax;
clear Valid391
Ratio(~Valid) = NaN;
CB = 1./Ratio;                                            % N2(C)/N2+(B), reported
if ReportEoN == 1
    EoN = polyval(EoNCoeff,Ratio);
    EoN(Ratio <= RatioRange(1) | Ratio >= RatioRange(2)) = NaN;
end

% Whole-electrode ratio per gate (sum over valid pixels, not a mean of ratios)
FrameRatio = nan(1,Nframes);
FrameEoN   = nan(1,Nframes);
for n = 1:Nframes
    v = Valid(:,:,n);
    if nnz(v) > 50
        a = A(:,:,n);  b = BFNS(:,:,n);
        FrameRatio(n) = (sum(b(v))/FilterFactor391)/(sum(a(v))/FilterFactor337);
        FrameEoN(n)   = polyval(EoNCoeff,FrameRatio(n));
    end
end
FrameCB = 1./FrameRatio;

RawFrameRatio = Tot391(:)'./Tot337(:)';
RawFrameRatio(~On) = NaN;

figure('Position',[100 100 900 700])
tiledlayout(3,1)
nexttile
plot(GateDelay,Tot337/1e6,'b','LineWidth',1.5); hold on
plot(GateDelay,Tot391/1e6,'r','LineWidth',1.5)
ylabel('Total counts / 10^6')
yyaxis right
plot(GateDelay,VGate/1e3,'k--'); ylabel('Voltage / kV')
legend('337 nm','391 nm','Applied voltage','Location','best')
xlabel('Gate delay / \mus')
title(sprintf('Emission and voltage at each gate (scope t = gate + %.1f \\mus)',ScopeTimeShift*1e6))
nexttile
plot(GateDelay,RawFrameRatio,'o-','LineWidth',1.5)
ylabel('Raw 391/337'); xlabel('Gate delay / \mus')
title('Electrode-integrated raw ratio (registered, no filter correction)')
nexttile
plot(GateDelay,FrameCB,'o-','LineWidth',1.5)
ylabel('N_2(C) / N_2^+(B)'); xlabel('Gate delay / \mus')
if ReportEoN == 1
    yyaxis right; plot(GateDelay,FrameEoN,'.-'); ylabel('E/N / Td')
end
title(sprintf('Electrode-integrated N_2(C)/N_2^+(B) emission ratio (%s, echelle-anchored)',ContaminationMode))
SaveFigure(gcf,OutDir.Ratio,'EmissionAndFrameRatio')
FrameTable = table(GateDelay(:),GateStart(:)*1e6,HalfWave(:),Tot337,Tot391,RawFrameRatio(:),FrameCB(:), ...
           'VariableNames',{'GateDelay_us','ScopeTime_us','HalfWave','Total337','Total391','Raw391over337','N2C_over_N2pB'});
if ReportEoN == 1, FrameTable.EoN_Td = FrameEoN(:); end
writetable(FrameTable,fullfile(OutDir.Ratio,'FrameEmissionRatio.csv'));

%% Section 7 - GIF of the ratio (or E/N) map and the current trace
if LogicGif == 1
    if isempty(GifFrames), GifFrames = 1:Nframes; end
    if exist(GifFile,'file'), delete(GifFile); end
    x_mm = ((1:Nx)-1)*ConversionFactor;
    y_mm = ((1:Ny)-1)*ConversionFactor;
    PlotWindow = [GateStart(1) GateStart(end)+GateWidth*1e-6];
    InWin = tScope >= PlotWindow(1) & tScope <= PlotWindow(2);
    h = figure('Color','w','Position',[100 100 640 680],'Name','BigHexEmissionRatio GIF (last frame)');
    ShowEoN = strcmpi(GifQuantity,'eon') && ReportEoN == 1;
    for nn = GifFrames
        clf(h)
        tiledlayout(h,3,1,'TileSpacing','compact','Padding','compact')
        nexttile([2 1])
        if ShowEoN
            Q = EoN(:,:,nn);  Lim = EoNLimits;  Cmap = parula;
            Lab = 'Reduced Electric Field / Td';  Ttl = 'EoN Data';
        else
            % high ratio = low field, so the map is flipped to keep bright = high field
            Q = CB(:,:,nn);   Lim = CBLimits;   Cmap = flipud(parula);
            Lab = 'N_2(C) / N_2^+(B) emission ratio';  Ttl = 'N_2(C) / N_2^+(B)';
        end
        imagesc(x_mm,y_mm,Q,'AlphaData',~isnan(Q)); axis image
        ax = gca; ax.Color = [0.12 0.12 0.3];
        colormap(Cmap); clim(Lim)
        if ~ShowEoN, ax.ColorScale = 'log'; end
        a = colorbar; a.Label.String = Lab;
        xlabel('x / mm'); ylabel('y / mm'); title(Ttl)
        text(x_mm(10),y_mm(end-15),sprintf('Gate = %d - %d \\mus',GateDelay(nn),GateDelay(nn)+GateWidth), ...
             'Color','w','FontSize',12)
        ax = gca; ax.FontSize = 12;
        nexttile
        plot(tScope(InWin)*1e3,IScope(InWin),'LineWidth',1.2); hold on
        xline(GateStart(nn)*1e3,'k'); xline((GateStart(nn)+GateWidth*1e-6)*1e3,'k')
        ylabel('Current / mA'); xlabel('Time / ms')
        xlim(PlotWindow*1e3)
        yyaxis right
        plot(tScope(InWin)*1e3,VScope(InWin)/1e3,'--','Color',[0.5 0.5 0.5])
        ylabel('Voltage / kV')
        ax = gca; ax.FontSize = 12;
        drawnow
        exportgraphics(h,GifFile,'Append',true);
    end
    fprintf('GIF written: %s (%d frames)\n',GifFile,numel(GifFrames));
end

%% Section 8 - Delayed ignition vs 391/337 ratio, per hexagon
% Superseded by Section 10 and off by default; still reports E/N, unanchored.
if LogicDelayed == 1
    Sk = load(SkelFile);
    Rings  = Sk.Rings;
    Ncells = numel(Rings);
    Band   = strel('disk',BandWidthPx);
    CellMask = cell(1,Ncells);
    RingBand = cell(1,Ncells);
    for k = 1:Ncells
        CellMask{k} = imdilate(imfill(Rings{k},'holes'),strel('disk',CellDilatePx));
        RingBand{k} = Rings{k};
    end

    % Per cell, per gate: lit share of the traced edge and summed signals
    LitShare = zeros(Ncells,Nframes);
    Cell337  = zeros(Ncells,Nframes);
    Cell391  = zeros(Ncells,Nframes);
    CellFNS  = zeros(Ncells,Nframes);
    for n = 1:Nframes
        if HalfWave(n) < 0
            Thr = IgnitionCountsNegative;
        else
            Thr = IgnitionCountsPositive;
        end
        Lit = bwareaopen(S337(:,:,n) > Thr,10);
        LitBand = imdilate(Lit,Band);
        a = double(S337(:,:,n));  b = double(S391(:,:,n));
        switch lower(ContaminationMode)
            case 'subtract', bf = b - Kappa*a;
            case 'fraction', bf = b*Factor391Percentage;
        end
        for k = 1:Ncells
            LitShare(k,n) = nnz(LitBand & RingBand{k})/nnz(RingBand{k});
            m = CellMask{k};
            Cell337(k,n) = sum(a(m));
            Cell391(k,n) = sum(b(m));
            CellFNS(k,n) = sum(bf(m));
        end
    end

    HWNames = {'Negative','Positive'};
    HWSign  = [-1 1];
    Results = table;
    for hw = 1:2
        g = find(HalfWave == HWSign(hw));
        % Ignition gate of each cell = first gate with LitFraction of its edge lit
        IgnIdx = nan(Ncells,1);
        for k = 1:Ncells
            j = find(LitShare(k,g) >= LitFraction,1,'first');
            if ~isempty(j), IgnIdx(k) = g(j); end
        end
        FirstIgn = min(IgnIdx);
        DelayGates = IgnIdx - FirstIgn;
        Delayed = DelayGates >= median(DelayGates,'omitnan') + DelayedFrames;
        AllLit  = max(IgnIdx);
        % Possible confounders: position in the field (vignetting, or a focus
        % difference between the two filters) and how bright the cell is
        Radius_mm = hypot(Sk.Centers(:,1) - mean(Sk.Centers(:,1)), ...
                          Sk.Centers(:,2) - mean(Sk.Centers(:,2)))*ConversionFactor;
        Signal337 = sum(Cell337(:,g),2);

        RatioHW    = nan(Ncells,1);   % whole half-wave
        RatioOnset = nan(Ncells,1);   % first OnsetFrames gates after own ignition
        RatioLate  = nan(Ncells,1);   % gates after every cell has ignited
        LateGates  = g(g >= AllLit);
        for k = 1:Ncells
            RatioHW(k) = CellRatio(CellFNS(k,g),Cell337(k,g),FilterFactor337,FilterFactor391);
            if ~isnan(IgnIdx(k))
                o = IgnIdx(k):min(IgnIdx(k)+OnsetFrames-1,g(end));
                RatioOnset(k) = CellRatio(CellFNS(k,o),Cell337(k,o),FilterFactor337,FilterFactor391);
            end
            if ~isempty(LateGates)
                RatioLate(k) = CellRatio(CellFNS(k,LateGates),Cell337(k,LateGates),FilterFactor337,FilterFactor391);
            end
        end
        T = table(repmat(HWNames(hw),Ncells,1),(1:Ncells)',Sk.Centers(:,1),Sk.Centers(:,2),Radius_mm,Signal337, ...
                  GateDelay(IgnIdx)',DelayGates*PictureDt,Delayed, ...
                  RatioHW,polyval(EoNCoeff,RatioHW),RatioOnset,polyval(EoNCoeff,RatioOnset), ...
                  RatioLate,polyval(EoNCoeff,RatioLate), ...
                  'VariableNames',{'HalfWave','Cell','Xc_px','Yc_px','Radius_mm','Signal337', ...
                                   'IgnitionGate_us','Delay_us','Delayed', ...
                                   'RatioHalfWave','EoNHalfWave_Td','RatioOnset','EoNOnset_Td','RatioLate','EoNLate_Td'});
        Results = [Results; T]; %#ok<AGROW>

        fprintf(['\n%s half-wave: first cells ignite at %d us, median %d us, all by %d us; ' ...
                 '%d of %d cells delayed (>= %d us after the median cell)\n'], ...
                HWNames{hw},GateDelay(FirstIgn),round(median(GateDelay(IgnIdx(~isnan(IgnIdx))))), ...
                GateDelay(AllLit),nnz(Delayed),Ncells,DelayedFrames*PictureDt);
        ok = ~isnan(T.Delay_us) & ~isnan(T.RatioHalfWave);
        fprintf('  Confounders: Spearman(delay,radius) = %.2f, (ratio,radius) = %.2f, (ratio,337 signal) = %.2f\n', ...
                corr(T.Delay_us(ok),Radius_mm(ok),'type','Spearman'), ...
                corr(T.RatioHalfWave(ok),Radius_mm(ok),'type','Spearman'), ...
                corr(T.RatioHalfWave(ok),Signal337(ok),'type','Spearman'));
        [rp,pp] = partialcorr(T.Delay_us(ok),T.RatioHalfWave(ok),[Radius_mm(ok) Signal337(ok)],'type','Spearman');
        fprintf('  Spearman(delay,half-wave ratio) controlling for radius and signal = %.2f (p = %.3g)\n',rp,pp);
        Metric = {'RatioHalfWave','RatioOnset','RatioLate'};
        for q = 1:3
            r = T.(Metric{q});
            e = r(~Delayed & ~isnan(r));  d = r(Delayed & ~isnan(r));
            if numel(d) >= 2 && numel(e) >= 2
                p   = ranksum(d,e);
                v   = ~isnan(r) & ~isnan(T.Delay_us);
                [rho,prho] = corr(T.Delay_us(v),r(v),'type','Spearman');
                fprintf(['  %-14s others %.4f +- %.4f (n=%d) | delayed %.4f +- %.4f (n=%d) | ranksum p = %.3g | ' ...
                         'Spearman(delay,ratio) = %.2f (p = %.3g)\n'], ...
                        Metric{q},mean(e),std(e),numel(e),mean(d),std(d),numel(d),p,rho,prho);
            else
                fprintf('  %-14s not enough delayed/early cells to compare\n',Metric{q});
            end
        end

        % Time-resolved ratio for the two groups
        GroupRatio = nan(2,numel(g));
        for j = 1:numel(g)
            for grp = 1:2
                sel = Delayed == (grp == 2);
                if any(sel)
                    GroupRatio(grp,j) = CellRatio(CellFNS(sel,g(j)),Cell337(sel,g(j)),FilterFactor337,FilterFactor391);
                end
            end
        end
        GroupSignal = [sum(Cell337(~Delayed,g),1); sum(Cell337(Delayed,g),1)];
        % Only plot the ratio where the group actually emits
        GroupRatio(GroupSignal < 0.05*max(GroupSignal,[],2)) = NaN;
        gLit = any(~isnan(GroupRatio),1);

        figure('Position',[100 100 1300 420])
        tiledlayout(1,3)
        nexttile
        CellDelayImg = nan(Ny,Nx);
        for k = 1:Ncells
            CellDelayImg(CellMask{k}) = T.Delay_us(k);
        end
        imagesc(CellDelayImg,'AlphaData',~isnan(CellDelayImg)); axis image off
        colormap(gca,turbo); a = colorbar; a.Label.String = 'Cell ignition delay / \mus';
        hold on
        for k = 1:Ncells
            text(Sk.Centers(k,1),Sk.Centers(k,2),num2str(k),'HorizontalAlignment','center', ...
                 'Color','k','FontWeight','bold','FontSize',8)
            if Delayed(k)
                [yy,xx] = find(Rings{k}); plot(xx,yy,'.m','MarkerSize',3)
            end
        end
        title(sprintf('%s half-wave (delayed cells outlined)',HWNames{hw}))
        nexttile
        hold on
        scatter(T.Delay_us,T.RatioHalfWave,40,'filled')
        scatter(T.Delay_us,T.RatioOnset,40,'^')
        scatter(T.Delay_us,T.RatioLate,40,'s')
        xlabel('Cell ignition delay / \mus'); ylabel('Corrected 391/337 ratio')
        legend('Whole half-wave','Own onset','After all lit','Location','best'); box on
        title('Ratio vs ignition delay, per hexagon')
        nexttile
        gp = g(gLit);
        yyaxis left
        plot(GateDelay(gp),GroupRatio(1,gLit),'b-o','MarkerSize',3); hold on
        plot(GateDelay(gp),GroupRatio(2,gLit),'r-o','MarkerSize',3)
        ylabel('Corrected 391/337 ratio')
        yyaxis right
        plot(GateDelay(gp),GroupSignal(1,gLit)/max(GroupSignal(:)),'b:')
        plot(GateDelay(gp),GroupSignal(2,gLit)/max(GroupSignal(:)),'r:')
        ylabel('337 signal (norm.)')
        xlabel('Gate delay / \mus'); legend('Other cells','Delayed cells','Location','best')
        title('Group ratio through the half-wave')
        SaveFigure(gcf,OutDir.Superseded,['DelayedIgnitionRatio_' HWNames{hw}])
    end
    writetable(Results,fullfile(OutDir.Superseded,'DelayedIgnitionRatio.csv'));
end

%% Section 9 - Delayed ignition vs 391/337 ratio, pixel by pixel
% Superseded by Section 10 and off by default; still reports E/N.
% Two thresholds: PixelOnCounts decides which pixels ever turn on in the
% half-wave, PixelDelayThreshold decides which of those count as delayed.
% Ratios are sums over pixels and gates of the smoothed, registered images
% (A, BFNS from Section 6), so dim pixels do not dominate through noise.
if LogicDelayedPixel == 1
    HWNames = {'Negative','Positive'};
    SkP = load(SkelFile);
    SkelAll = SkP.Skeleton;     % traced hexagon edges + outer perimeter
    [XX,YY] = meshgrid(1:Nx,1:Ny);
    Ctr = [mean(XX(any(A > PixelOnCounts(2),3))) mean(YY(any(A > PixelOnCounts(2),3)))];
    RadMap = hypot(XX - Ctr(1),YY - Ctr(2))*ConversionFactor;   % mm from electrode centre
    PixSummary = table;
    for h = 1:numel(HumpStart)
        g  = HumpStart(h):HumpEnd(h);
        hw = 1 + (HalfWave(g(1)) > 0);
        Ah = A(:,:,g);  Fh = BFNS(:,:,g);
        PeakA = max(Ah,[],3);
        EverOn = PeakA >= PixelOnCounts(hw);

        % Ignition gate per pixel, relative to its own maximum
        Above = Ah >= PixelIgnFraction*PeakA;
        [~,IgnLocal] = max(Above,[],3);               % first gate above threshold
        IgnTime = GateDelay(g(1)) + (IgnLocal - 1)*PictureDt;
        IgnTime(~EverOn) = NaN;
        % Delay is only judged on the electrode edge, where the discharge ignites.
        % Pixels further out light later because the plasma expands, so they
        % take the label of their nearest edge pixel instead of their own time.
        Edge    = EverOn & imdilate(SkelAll,strel('disk',BandWidthPx));
        RefTime = prctile(IgnTime(Edge),PixelRefPercentile);
        Delay   = IgnTime - RefTime;
        EdgeDelayed = Edge & Delay >= PixelDelayThreshold;
        [~,Nearest] = bwdist(Edge);
        Delayed = EverOn & EdgeDelayed(Nearest);
        Prompt  = EverOn & ~Delayed;

        % Per-pixel ratios: whole half-wave, and first OnsetFrames gates after own ignition
        RatioHWPix = (sum(Fh,3)/FilterFactor391)./(sum(Ah,3)/FilterFactor337);
        k3 = reshape(1:numel(g),1,1,[]);
        OnWin = k3 >= IgnLocal & k3 < IgnLocal + OnsetFrames;
        RatioOnPix = (sum(Fh.*OnWin,3)/FilterFactor391)./(sum(Ah.*OnWin,3)/FilterFactor337);
        RatioHWPix(~EverOn) = NaN;  RatioOnPix(~EverOn) = NaN;

        GroupR = @(M,F,Am) (sum(F(M))/FilterFactor391)/(sum(Am(M))/FilterFactor337);
        SA = sum(Ah,3); SF = sum(Fh,3); OA = sum(Ah.*OnWin,3); OF = sum(Fh.*OnWin,3);
        R = [GroupR(Prompt,SF,SA) GroupR(Delayed,SF,SA); GroupR(Prompt,OF,OA) GroupR(Delayed,OF,OA)];

        % Position-controlled: delayed minus prompt within radial strata, signal weighted
        Edges = linspace(0,max(RadMap(EverOn)),RadiusBins+1);
        dHW = nan(1,RadiusBins); wts = zeros(1,RadiusBins);
        for b = 1:RadiusBins
            inb = RadMap >= Edges(b) & RadMap < Edges(b+1);
            if nnz(inb & Delayed) > 200 && nnz(inb & Prompt) > 200
                dHW(b) = GroupR(inb & Delayed,SF,SA)/GroupR(inb & Prompt,SF,SA) - 1;
                wts(b) = sum(SA(inb & EverOn));
            end
        end
        dStrat = sum(dHW(wts>0).*wts(wts>0))/sum(wts);

        fprintf(['\n%s half-wave (pixels): %d turn on, %d of %d edge pixels delayed (>= %g us after the ' ...
                 '%dth-percentile edge ignition at %.0f us), %d pixels labelled delayed\n'], ...
                 HWNames{hw},nnz(EverOn),nnz(EdgeDelayed),nnz(Edge),PixelDelayThreshold, ...
                 PixelRefPercentile,RefTime,nnz(Delayed));
        Re = [GroupR(Edge & ~EdgeDelayed,SF,SA) GroupR(EdgeDelayed,SF,SA); ...
              GroupR(Edge & ~EdgeDelayed,OF,OA) GroupR(EdgeDelayed,OF,OA)];
        fprintf('  Edge pixels only     : half-wave %+.1f%%, own onset %+.1f%% (delayed vs prompt)\n', ...
                100*(Re(1,2)/Re(1,1)-1),100*(Re(2,2)/Re(2,1)-1));
        fprintf('  Whole half-wave ratio: prompt %.4f (%.0f Td) | delayed %.4f (%.0f Td) | %+.1f%%\n', ...
                R(1,1),polyval(EoNCoeff,R(1,1)),R(1,2),polyval(EoNCoeff,R(1,2)),100*(R(1,2)/R(1,1)-1));
        fprintf('  Own-onset ratio      : prompt %.4f (%.0f Td) | delayed %.4f (%.0f Td) | %+.1f%%\n', ...
                R(2,1),polyval(EoNCoeff,R(2,1)),R(2,2),polyval(EoNCoeff,R(2,2)),100*(R(2,2)/R(2,1)-1));
        fprintf('  Whole half-wave, same-radius comparison: %+.1f%% (per stratum: %s)\n', ...
                100*dStrat,sprintf('%+.1f ',100*dHW));
        PixSummary = [PixSummary; table(HWNames(hw),nnz(EverOn),nnz(Delayed),RefTime, ...
            R(1,1),R(1,2),R(2,1),R(2,2),polyval(EoNCoeff,R(1,1)),polyval(EoNCoeff,R(1,2)), ...
            polyval(EoNCoeff,R(2,1)),polyval(EoNCoeff,R(2,2)),dStrat, ...
            'VariableNames',{'HalfWave','PixelsOn','PixelsDelayed','RefIgnition_us', ...
            'RatioPrompt','RatioDelayed','OnsetRatioPrompt','OnsetRatioDelayed', ...
            'EoNPrompt_Td','EoNDelayed_Td','OnsetEoNPrompt_Td','OnsetEoNDelayed_Td','SameRadiusDiff'})]; %#ok<AGROW>

        % Group ratio through the half-wave
        GR = nan(2,numel(g));  GS = zeros(2,numel(g));
        for j = 1:numel(g)
            a = Ah(:,:,j); f = Fh(:,:,j);
            GS(:,j) = [sum(a(Prompt)); sum(a(Delayed))];
            GR(:,j) = [GroupR(Prompt,f,a); GroupR(Delayed,f,a)];
        end
        GR(GS < 0.05*max(GS,[],2)) = NaN;

        figure('Position',[50 50 1500 820])
        tiledlayout(2,3)
        nexttile
        imagesc(Delay,'AlphaData',EverOn); axis image off
        colormap(gca,turbo); c = colorbar; c.Label.String = 'Ignition delay / \mus';
        title(sprintf('%s: ignition delay (ref %.0f \\mus)',HWNames{hw},RefTime))
        nexttile
        Cls = zeros(Ny,Nx); Cls(Prompt) = 1; Cls(Delayed) = 2;
        imagesc(Cls); axis image off; colormap(gca,[1 1 1; 0.3 0.5 0.9; 0.9 0.2 0.2])
        title(sprintf('Blue prompt, red delayed (>= %g \\mus)',PixelDelayThreshold))
        nexttile
        E = polyval(EoNCoeff,RatioHWPix);
        imagesc(E,'AlphaData',EverOn & ~isnan(E)); axis image off
        colormap(gca,parula); clim(EoNLimits); c = colorbar; c.Label.String = 'Half-wave E/N / Td';
        title('Half-wave integrated E/N')
        nexttile
        Bins = linspace(EoNLimits(1),EoNLimits(2),60);
        histogram(E(Prompt),Bins,'Normalization','probability','FaceColor',[0.3 0.5 0.9]); hold on
        histogram(E(Delayed),Bins,'Normalization','probability','FaceColor',[0.9 0.2 0.2])
        xlabel('Half-wave E/N / Td'); ylabel('Pixel fraction'); legend('Prompt','Delayed')
        title('Per-pixel E/N distribution')
        nexttile
        scatter(Delay(EverOn),RatioOnPix(EverOn),2,RadMap(EverOn),'filled','MarkerFaceAlpha',0.2)
        c = colorbar; c.Label.String = 'Radius / mm'; xline(PixelDelayThreshold,'k--')
        xlabel('Ignition delay / \mus'); ylabel('Own-onset ratio'); ylim([0 2*R(2,1)])
        title('Onset ratio vs delay, per pixel')
        nexttile
        yyaxis left
        plot(GateDelay(g),GR(1,:),'b-o','MarkerSize',3); hold on
        plot(GateDelay(g),GR(2,:),'r-o','MarkerSize',3); ylabel('Corrected 391/337')
        yyaxis right
        plot(GateDelay(g),GS(1,:)/max(GS(:)),'b:'); plot(GateDelay(g),GS(2,:)/max(GS(:)),'r:')
        ylabel('337 signal (norm.)'); xlabel('Gate delay / \mus')
        legend('Prompt','Delayed','Location','best'); title('Group ratio through the half-wave')
        SaveFigure(gcf,OutDir.Superseded,['DelayedPixelRatio_' HWNames{hw}])
    end
    writetable(PixSummary,fullfile(OutDir.Superseded,'DelayedPixelRatio.csv'));
end

%% Section 10 - Per-side rectangles: ignition delay vs E/N, and the open perimeter ratio
% Every stamped hexagon (PerimeterSkeleton.mat) is split into its 6 sides. A
% side owns the 60 degree wedge of its hexagon between the centre-to-vertex
% lines that bound it, so neighbouring sides never share pixels. Within that
% wedge the side's rectangle runs SideOutMM over the foil edge and SideDepthMM
% into the opening, which holds the whole expansion region of that side (the
% Fig 3 extent rectangle). Corners (SideTrim of the length at each end) are
% dropped because the plasma from two sides merges there.
%
% Ignition: the same lit test as PerimeterCodev2.m (337 above the half-wave
% threshold, despeckled, dilated by BandWidthPx), applied to a thin line along
% the side. A side ignites at the first gate with SideLitFraction of its line lit.
%
% E/N per side is compared in three windows:
%   Matched - gates after SideMatchedPct % of sides have ignited, i.e. the same
%             gates for every side. The ratio drifts through each half-wave, so
%             this is the fair comparison.
%   Onset   - the first OnsetFrames gates after the side's own ignition
%   HalfWave- the whole half-wave
if LogicSides == 1
    SkS = load(SkelFile);
    Centers = SkS.Centers;
    Ncells  = size(Centers,1);
    Rpx     = SkS.SideMM/ConversionFactor;      % circumradius = side length (stamp geometry)
    Apothem = Rpx*cosd(30);
    InPx    = SideDepthMM/ConversionFactor;
    OutPx   = SideOutMM/ConversionFactor;
    if InPx >= Apothem
        error('SideDepthMM reaches past the hexagon centre (apothem %.2f mm).',Apothem*ConversionFactor);
    end
    Nsides = 6*Ncells;
    [XX,YY] = meshgrid(1:Nx,1:Ny);
    Phis    = 60*(0:5)' + 30;                     % direction of each side's midpoint

    % Refit each hexagon: the time-summed 337 light peaks at the electrode edge
    % (Fig 10C), so find that peak across each side and solve for the shift of
    % the centre (dx, dy) and the change in apothem that best explains all six:
    %   peak depth of side j = dA - dx*cos(phi_j) - dy*sin(phi_j)
    ApothemK = Apothem*ones(Ncells,1);
    if RefineHexagons == 1
        LightSum = double(Sum337);
        dgrid = -12:20;
        Shift = zeros(Ncells,3);
        for k = 1:Ncells
            rx = XX - Centers(k,1);  ry = YY - Centers(k,2);
            th = mod(atan2d(ry,rx),360);
            Pk = nan(6,1);
            for j = 1:6
                Wedge = th >= 60*(j-1) & th < 60*j;
                d = round(Apothem - (rx*cosd(Phis(j)) + ry*sind(Phis(j))));
                s = -rx*sind(Phis(j)) + ry*cosd(Phis(j));
                m = Wedge & abs(s) <= (0.5 - SideTrim)*Rpx & d >= dgrid(1) & d <= dgrid(end);
                p = accumarray(d(m) - dgrid(1) + 1,LightSum(m),[numel(dgrid) 1],@mean);
                [~,i] = max(p);  Pk(j) = dgrid(i);
            end
            Shift(k,:) = ([ones(6,1) -cosd(Phis) -sind(Phis)] \ Pk)';
            if any(abs(Shift(k,:)) > MaxRefinePx)
                warning('Hexagon %d refit rejected (%.1f, %.1f, %.1f px); keeping the stamp.',k,Shift(k,:));
                Shift(k,:) = 0;
            end
        end
        Centers  = Centers + Shift(:,2:3);
        ApothemK = Apothem - Shift(:,1);
        fprintf(['\nHexagon refit: mean centre shift dx = %.1f px, dy = %.1f px; mean side ' ...
                 'length %.2f mm (stamp %.2f mm)\n'],mean(Shift(:,2)),mean(Shift(:,3)), ...
                 mean(ApothemK)/cosd(30)*ConversionFactor,SkS.SideMM);
    end
    % Two geometries: 1 = the stamps as saved (what PerimeterCodev2.m and the
    % paper's Fig 13 used), 2 = the refitted hexagons. Geometry 2 drives the
    % per-side E/N statistics; both are kept so Fig13Compare.m can show old vs new.
    StampApothem = Apothem*ones(Ncells,1);
    G1 = BuildSideGeometry(SkS.Centers,StampApothem,Ny,Nx,SideTrim,InPx,OutPx,SideEdgeHalfPx);
    G1.Skel = SkS.Skeleton;
    G2 = BuildSideGeometry(Centers,ApothemK,Ny,Nx,SideTrim,InPx,OutPx,SideEdgeHalfPx);
    G2.Skel = G2.HexRings;                       % refitted hexagons + the traced outer arcs
    for m = 1:numel(SkS.FreeRings), G2.Skel = G2.Skel | SkS.FreeRings{m}; end
    if RefineHexagons ~= 1, G2 = G1; end
    Geo = [G1 G2];
    RpxK      = ApothemK/cosd(30);
    RectLabel = Geo(2).RectLabel;   EdgeLabel = Geo(2).EdgeLabel;
    SideMid   = Geo(2).SideMid;     SkelAll   = Geo(2).Skel;
    SideCell  = repelem((1:Ncells)',6);
    SideNum   = repmat((1:6)',Ncells,1);
    DepthMM   = Geo(2).DepthPx*ConversionFactor;   % depth of each extent bin from the edge
    OuterArcMM = nnz(any(cat(3,SkS.FreeRings{:}),3))*ConversionFactor;   % traced outer perimeter (Pout, Eq 14)
    save(fullfile(OutDir.Cache,'HexagonGeometryRefined.mat'),'Centers','ApothemK','RpxK','SkelAll', ...
         'RectLabel','EdgeLabel','SideCell','SideNum','SideMid','ConversionFactor');

    % QC: rectangles over the time-summed 337 image
    figure('Position',[50 50 900 800])
    B = labeloverlay(rescale(Sum337),RectLabel,'Colormap',repmat(lines(6),Ncells,1),'Transparency',0.55);
    image(B); axis image off; hold on
    for id = 1:Nsides
        text(SideMid(id,1),SideMid(id,2),num2str(SideNum(id)),'Color','y','FontSize',6, ...
             'HorizontalAlignment','center')
    end
    title('Side rectangles (numbers = side within its hexagon)')
    SaveFigure(gcf,OutDir.Sides,'SideRectangles_QC')

    % Per gate, for both geometries: lit share of each edge line, summed signals
    % in each rectangle, mean 337 intensity at each depth (the Fig 3 profile),
    % and the open perimeter ratio. Plasma area and total intensity use the same
    % threshold mask as PerimeterCodev2.m and do not depend on the geometry.
    Nd = numel(Geo(2).DepthPx);
    SideLit  = zeros(Nsides,Nframes,2);
    SideProf = zeros(Nsides,Nd,Nframes,2);
    SideProfFNS = zeros(Nsides,Nd,Nframes);   % refitted geometry only
    MatchStart  = nan(1,numel(HumpStart));    % first gate of each half-wave's matched window
    Zeta     = zeros(2,Nframes);
    Side337  = zeros(Nsides,Nframes);
    SideFNS  = zeros(Nsides,Nframes);
    PlasmaArea = zeros(1,Nframes);   % mm^2
    ArcMask    = any(cat(3,SkS.FreeRings{:}),3);   % traced outer perimeter
    ArcLit     = zeros(1,Nframes);   % lit share of the outer perimeter
    TotalInt   = zeros(1,Nframes);   % 337 counts inside the plasma mask
    for n = 1:Nframes
        if HalfWave(n) < 0, Thr = IgnitionCountsNegative; else, Thr = IgnitionCountsPositive; end
        Mask    = bwareaopen(S337(:,:,n) > Thr,10);
        LitBand = imdilate(Mask,strel('disk',BandWidthPx));
        a = double(S337(:,:,n));  b = AnchorGain*double(S391(:,:,n));
        switch lower(ContaminationMode)
            case 'subtract', bf = b - Kappa*a;
            case 'fraction', bf = b*Factor391Percentage;
        end
        ArcLit(n)     = nnz(LitBand & ArcMask)/nnz(ArcMask);
        PlasmaArea(n) = nnz(Mask)*ConversionFactor^2;
        TotalInt(n)   = sum(a(Mask));
        for q = 1:2
            G = Geo(q);
            Zeta(q,n) = nnz(LitBand & G.Skel)/nnz(G.Skel);
            em = G.EdgeLabel > 0;
            SideLit(:,n,q) = accumarray(G.EdgeLabel(em),LitBand(em),[Nsides 1])./G.EdgeN;
            dm = G.DepthLabel > 0;
            SideProf(:,:,n,q) = accumarray([G.RectLabel(dm) G.DepthLabel(dm)],a(dm),[Nsides Nd]) ...
                                ./max(G.DepthN,1);
            if q == 2
                % same profile for the corrected N2+ signal, so the ratio (and
                % E/N) can be compared at the same distance from the edge
                SideProfFNS(:,:,n) = accumarray([G.RectLabel(dm) G.DepthLabel(dm)],bf(dm),[Nsides Nd]) ...
                                     ./max(G.DepthN,1);
            end
        end
        rm = RectLabel > 0;
        Side337(:,n) = accumarray(RectLabel(rm),a(rm),[Nsides 1]);
        SideFNS(:,n) = accumarray(RectLabel(rm),bf(rm),[Nsides 1]);
    end

    % Plasma extent per side and gate (Fig 3 method): walking in from the edge,
    % the first depth where the mean 337 intensity falls below ExtentCounts.
    SideExtent = nan(Nsides,Nframes,2);
    for q = 1:2
        for n = 1:Nframes
            if HalfWave(n) < 0, Thr = ExtentCounts(1); else, Thr = ExtentCounts(2); end
            P  = SideProf(:,:,n,q);
            In = P(:,DepthMM >= 0) >= Thr;
            Dm = DepthMM(DepthMM >= 0);
            for id = 1:Nsides
                if In(id,1)
                    j = find(~In(id,:),1,'first');
                    if isempty(j), SideExtent(id,n,q) = Dm(end); else, SideExtent(id,n,q) = Dm(j); end
                else
                    SideExtent(id,n,q) = 0;
                end
            end
        end
    end
    SideLitOld = SideLit(:,:,1);
    SideLit    = SideLit(:,:,2);
    ZetaOld    = Zeta(1,:);
    Zeta       = Zeta(2,:);
    SideR  = @(F,A) (sum(F,2)/FilterFactor391)./(sum(A,2)/FilterFactor337);   % N2+(B)/N2(C)
    SideCB = @(F,A) 1./SideR(F,A);                                             % N2(C)/N2+(B)

    HWNames = {'Negative','Positive'};
    SideTable = table;
    DelayBoth = nan(Nsides,2);
    DelayedBoth = false(Nsides,2);
    for h = 1:numel(HumpStart)
        g  = HumpStart(h):HumpEnd(h);
        hw = 1 + (HalfWave(g(1)) > 0);
        IgnIdx = nan(Nsides,1);
        for id = 1:Nsides
            j = find(SideLit(id,g) >= SideLitFraction,1,'first');
            if ~isempty(j), IgnIdx(id) = g(j); end
        end
        Ignited  = ~isnan(IgnIdx);
        IgnTime  = nan(Nsides,1);  IgnTime(Ignited) = GateDelay(IgnIdx(Ignited));
        IgnVolt  = nan(Nsides,1);  IgnVolt(Ignited) = VGate(IgnIdx(Ignited));
        Delay    = IgnTime - median(IgnTime,'omitnan');
        Delayed  = Ignited & Delay >= SideDelayThreshold;
        MatchIdx = round(prctile(IgnIdx(Ignited),SideMatchedPct));
        MatchStart(h) = MatchIdx;
        gm = MatchIdx:g(end);

        CBMatched = SideCB(SideFNS(:,gm),Side337(:,gm));
        CBHW      = SideCB(SideFNS(:,g),Side337(:,g));
        CBOnset   = nan(Nsides,1);
        for id = find(Ignited)'
            o = IgnIdx(id):min(IgnIdx(id)+OnsetFrames-1,g(end));
            CBOnset(id) = SideCB(SideFNS(id,o),Side337(id,o));
        end
        Signal = sum(Side337(:,g),2);
        Radius = hypot(SideMid(:,1)-mean(Centers(:,1)),SideMid(:,2)-mean(Centers(:,2)))*ConversionFactor;
        DelayBoth(:,hw) = Delay;  DelayedBoth(:,hw) = Delayed;

        T = table(repmat(HWNames(hw),Nsides,1),SideCell,SideNum,SideMid(:,1),SideMid(:,2),Radius,Signal, ...
                  IgnTime,IgnVolt,Delay,Delayed,CBMatched,CBOnset,CBHW, ...
                  'VariableNames',{'HalfWave','Cell','Side','Xmid_px','Ymid_px','Radius_mm','Signal337', ...
                  'IgnitionTime_us','IgnitionVoltage_V','Delay_us','Delayed', ...
                  'CBMatched','CBOnset','CBHalfWave'});
        if ReportEoN == 1
            T.EoNMatched_Td  = polyval(EoNCoeff,1./CBMatched);
            T.EoNOnset_Td    = polyval(EoNCoeff,1./CBOnset);
            T.EoNHalfWave_Td = polyval(EoNCoeff,1./CBHW);
        end
        SideTable = [SideTable; T]; %#ok<AGROW>

        fprintf(['\n%s half-wave, %d sides: %d ignite (median %.0f us, %.0f V); %d delayed (>= %g us ' ...
                 'after the median); matched window %d-%d us\n'],HWNames{hw},Nsides,nnz(Ignited), ...
                 median(IgnTime,'omitnan'),median(IgnVolt,'omitnan'),nnz(Delayed),SideDelayThreshold, ...
                 GateDelay(gm(1)),GateDelay(gm(end)));
        % A geometry error shows up as delays piling onto one side number
        fprintf('  Delayed sides by side number 1-6: %s\n',sprintf('%d ',accumarray(SideNum,Delayed)));
        Metric = {'CBMatched','CBOnset','CBHalfWave'};
        for q = 1:3
            r = T.(Metric{q});
            v = Ignited & ~isnan(r);
            e = r(v & ~Delayed);  d = r(v & Delayed);
            [rho,prho] = corr(Delay(v),r(v),'type','Spearman');
            if numel(d) >= 3
                p = ranksum(d,e);
                fprintf(['  %-10s N2(C)/N2+(B): others %.2f (n=%d) | delayed %.2f (n=%d) | %+.1f%% | ' ...
                         'ranksum p = %.3g | Spearman(delay,ratio) %.2f (p = %.3g)\n'],Metric{q}, ...
                         median(e),numel(e),median(d),numel(d),100*(median(d)/median(e)-1),p,rho,prho);
            end
        end
        % Within-hexagon: each delayed side against the other sides of its own
        % hexagon (removes position, focus and cell-to-cell differences)
        Diff = [];
        for k = unique(SideCell(Delayed))'
            in = SideCell == k & Ignited;
            for id = find(in & Delayed)'
                Others = in & ~Delayed;
                if any(Others)
                    Diff(end+1) = CBMatched(id)/median(CBMatched(Others)) - 1; %#ok<AGROW>
                end
            end
        end
        if numel(Diff) >= 3
            fprintf('  Within-hexagon, matched window: delayed side vs rest of its hexagon %+.1f%% (median, n = %d, signrank p = %.3g)\n', ...
                    100*median(Diff),numel(Diff),signrank(Diff));
        end
        v = Ignited & ~isnan(CBMatched);
        [rp,pp] = partialcorr(Delay(v),CBMatched(v),[Radius(v) Signal(v)],'type','Spearman');
        fprintf('  Spearman(delay, matched ratio) controlling for radius and signal: %.2f (p = %.3g)\n',rp,pp);

        % Figures: side maps of delay and matched ratio, and the scatter
        figure('Position',[50 50 1500 480])
        tiledlayout(1,3)
        nexttile
        M = nan(Ny,Nx);  M(rm) = Delay(RectLabel(rm));
        imagesc(M,'AlphaData',~isnan(M)); axis image off
        colormap(gca,turbo); c = colorbar; c.Label.String = 'Side ignition delay / \mus';
        hold on
        for id = find(Delayed)'
            plot(SideMid(id,1),SideMid(id,2),'mo','MarkerSize',6,'LineWidth',1.5)
        end
        title(sprintf('%s: side ignition delay (circles = delayed)',HWNames{hw}))
        nexttile
        M = nan(Ny,Nx);  M(rm) = CBMatched(RectLabel(rm));
        imagesc(M,'AlphaData',~isnan(M)); axis image off
        colormap(gca,flipud(parula)); c = colorbar; c.Label.String = 'N_2(C)/N_2^+(B), matched window';
        hold on
        for id = find(Delayed)'
            plot(SideMid(id,1),SideMid(id,2),'mo','MarkerSize',6,'LineWidth',1.5)
        end
        title(sprintf('N_2(C)/N_2^+(B) per side, gates %d-%d \\mus',GateDelay(gm(1)),GateDelay(gm(end))))
        nexttile
        hold on
        scatter(Delay(~Delayed),CBMatched(~Delayed),25,[0.3 0.5 0.9],'filled')
        scatter(Delay(Delayed),CBMatched(Delayed),25,[0.9 0.2 0.2],'filled')
        scatter(Delay(Ignited),CBOnset(Ignited),12,[0.5 0.5 0.5])
        xline(SideDelayThreshold,'k--')
        xlabel('Side ignition delay / \mus'); ylabel('N_2(C) / N_2^+(B)'); box on
        legend('Matched, other sides','Matched, delayed sides','Own onset (all sides)','Location','best')
        title('N_2(C)/N_2^+(B) vs ignition delay, per side')
        SaveFigure(gcf,OutDir.Sides,['SideIgnitionRatio_' HWNames{hw}])
    end
    writetable(SideTable,fullfile(OutDir.Sides,'SideIgnitionRatio.csv'));

    % Are the same sides late in both polarities?
    v = all(~isnan(DelayBoth),2);
    [rs,ps] = corr(DelayBoth(v,1),DelayBoth(v,2),'type','Spearman');
    Tab = [nnz(DelayedBoth(v,1) & DelayedBoth(v,2)) nnz(DelayedBoth(v,1) & ~DelayedBoth(v,2)); ...
           nnz(~DelayedBoth(v,1) & DelayedBoth(v,2)) nnz(~DelayedBoth(v,1) & ~DelayedBoth(v,2))];
    [~,pf] = fishertest(Tab);
    fprintf(['\nSame sides late in both half-waves? Spearman(delay neg, delay pos) = %.2f (p = %.3g, n = %d); ' ...
             'delayed in both: %d, neg only: %d, pos only: %d (Fisher p = %.3g)\n'], ...
             rs,ps,nnz(v),Tab(1,1),Tab(1,2),Tab(2,1),pf);
    figure
    scatter(DelayBoth(v,1),DelayBoth(v,2),25,'filled','MarkerFaceAlpha',0.5); box on
    xlabel('Negative half-wave delay / \mus'); ylabel('Positive half-wave delay / \mus')
    title(sprintf('Per-side ignition delay, Spearman %.2f (p = %.2g)',rs,ps))
    SaveFigure(gcf,OutDir.Sides,'SideDelayNegVsPos')

    % Open perimeter ratio through the period and against the applied voltage
    figure('Position',[50 50 1100 420])
    tiledlayout(1,2)
    nexttile
    plot(GateDelay,Zeta,'LineWidth',1.5); hold on
    yyaxis right; plot(GateDelay,VGate/1e3,'k--'); ylabel('Voltage / kV')
    yyaxis left; ylabel('Open perimeter ratio \zeta'); xlabel('Gate delay / \mus')
    nexttile; hold on
    for h = 1:numel(HumpStart)
        g = HumpStart(h):HumpEnd(h);
        plot(VGate(g),Zeta(g),'o-','MarkerSize',3)
    end
    xlabel('U / V'); ylabel('\zeta'); box on
    legend(HWNames(1 + (HalfWave(HumpStart) > 0)),'Location','best')
    title(sprintf('\\zeta vs applied voltage (scope t = gate + %.0f \\mus)',ScopeTimeShift*1e6))
    SaveFigure(gcf,OutDir.Perimeter,'OpenPerimeterRatio')
    writetable(table(GateDelay(:),VGate(:),HalfWave(:),Zeta(:), ...
               'VariableNames',{'GateDelay_us','Voltage_V','HalfWave','Zeta'}), ...
               fullfile(OutDir.Perimeter,'OpenPerimeterRatio.csv'));

    % Cache everything derived from the images, so comparisons and figures
    % (Fig13Compare.m) can be redone in seconds without touching the stacks.
    % Index 1 / "Old" = stamped hexagons, 2 / no suffix = refitted hexagons.
    Derived = struct('GateDelay',GateDelay,'GateWidth',GateWidth,'HalfWave',HalfWave, ...
        'HumpStart',HumpStart,'HumpEnd',HumpEnd,'ScopeTimeShift',ScopeTimeShift, ...
        'OldTriggerOffset',1.71628e-6,'ScopeFile',ScopeFile,'Fdrive',Fdrive, ...
        'AmpGain',AmpGain,'CurrentGain',CurrentGain,'ConversionFactor',ConversionFactor, ...
        'Zeta',Zeta,'ZetaOld',ZetaOld,'PlasmaArea',PlasmaArea,'TotalInt',TotalInt, ...
        'Tot337',Tot337(:)','Tot391',Tot391(:)','FrameRatio',FrameRatio,'FrameCB',FrameCB, ...
        'AnchorGain',AnchorGain,'EchelleRaw',EchelleRaw,'FNSFraction',Factor391Percentage,'Kappa',Kappa, ...
        'ArcLit',ArcLit,'ArcLenMM',nnz(bwskel(ArcMask))*ConversionFactor, ...
        'SideLenMM',repelem(RpxK*ConversionFactor,6),'ReportEoN',ReportEoN, ...
        'SideLit',SideLit,'SideLitOld',SideLitOld,'Side337',Side337,'SideFNS',SideFNS, ...
        'SideProf',SideProf,'SideProfFNS',SideProfFNS,'MatchStart',MatchStart, ...
        'SideExtent',SideExtent,'DepthMM',DepthMM,'ExtentCounts',ExtentCounts, ...
        'Ncells',Ncells,'SideLengthMM',mean(RpxK)*ConversionFactor, ...
        'SideLengthOldMM',SkS.SideMM,'OuterArcMM',OuterArcMM, ...
        'SideCell',SideCell,'SideNum',SideNum,'SideMid',SideMid,'SideTable',SideTable, ...
        'FilterFactor337',FilterFactor337,'FilterFactor391',FilterFactor391, ...
        'ContaminationMode',ContaminationMode,'EoNCoeff',EoNCoeff, ...
        'IgnitionCounts',[IgnitionCountsNegative IgnitionCountsPositive]);
    save(DerivedFile,'-struct','Derived','-v7.3');
    fprintf('Derived image data cached: %s\n',DerivedFile);
end

%% Local functions
function SaveFigure(fig,Folder,Name)
% Save a figure as .png and .fig and leave it open, titled with its file name.
fig.Name = Name;
fig.Visible = 'on';
exportgraphics(fig,fullfile(Folder,[Name '.png']))
savefig(fig,fullfile(Folder,[Name '.fig']))
end

function G = BuildSideGeometry(Centers,ApothemK,Ny,Nx,SideTrim,InPx,OutPx,EdgeHalfPx)
% Label maps for the 6 side rectangles of every hexagon (see Section 10).
%   RectLabel  - side index of each rectangle pixel
%   EdgeLabel  - side index of each pixel on the thin edge line
%   DepthLabel - depth bin (1 px bins from -OutPx to InPx) of each rectangle pixel
Ncells = size(Centers,1);
Nsides = 6*Ncells;
Phis = 60*(0:5)' + 30;
[XX,YY] = meshgrid(1:Nx,1:Ny);
G.DepthPx    = -floor(OutPx):floor(InPx);
G.RectLabel  = zeros(Ny,Nx);
G.EdgeLabel  = zeros(Ny,Nx);
G.DepthLabel = zeros(Ny,Nx);
G.SideMid    = zeros(Nsides,2);
G.HexRings   = false(Ny,Nx);
for k = 1:Ncells
    Rk = ApothemK(k)/cosd(30);
    rx = XX - Centers(k,1);
    ry = YY - Centers(k,2);
    Near = hypot(rx,ry) <= Rk + OutPx + 2;
    th = mod(atan2d(ry,rx),360);              % image coordinates, y down, as stamped
    Ang = (0:5)*60;
    G.HexRings = G.HexRings | bwperim(poly2mask(Centers(k,1) + Rk*cosd(Ang),Centers(k,2) + Rk*sind(Ang),Ny,Nx));
    for j = 1:6
        Wedge = th >= 60*(j-1) & th < 60*j;
        d = ApothemK(k) - (rx*cosd(Phis(j)) + ry*sind(Phis(j)));   % depth inward from the edge
        s = -rx*sind(Phis(j)) + ry*cosd(Phis(j));                   % position along the side
        Along = abs(s) <= (0.5 - SideTrim)*Rk;
        id = 6*(k-1) + j;
        In = Near & Wedge & Along & round(d) >= G.DepthPx(1) & round(d) <= G.DepthPx(end);
        G.RectLabel(In)  = id;
        G.DepthLabel(In) = round(d(In)) - G.DepthPx(1) + 1;
        G.EdgeLabel(Near & Wedge & Along & abs(d) <= EdgeHalfPx) = id;
        G.SideMid(id,:) = Centers(k,:) + ApothemK(k)*[cosd(Phis(j)) sind(Phis(j))];
    end
end
G.EdgeN  = accumarray(G.EdgeLabel(G.EdgeLabel>0),1,[Nsides 1]);
dm = G.DepthLabel > 0;
G.DepthN = accumarray([G.RectLabel(dm) G.DepthLabel(dm)],1,[Nsides numel(G.DepthPx)]);
end

function Img = ReadPimaxCSV(f)
% LightField CSV export: header 'Y,X,I', then one row per pixel, X fastest.
fid = fopen(f,'r');
C = textscan(fid,'%f%f%f','Delimiter',',','HeaderLines',1);
fclose(fid);
v = C{3};
if numel(v) < 1024^2
    warning('%s has %d of %d pixels; padding with zeros.',f,numel(v),1024^2);
    v(end+1:1024^2) = 0;
end
Img = reshape(v(1:1024^2),1024,1024)';
end

function R = CellRatio(FNS,S337,F337,F391)
% Ratio of summed signals (not a mean of per-gate ratios), filter corrected
s = sum(S337(:));
if s <= 0
    R = NaN;
else
    R = (sum(FNS(:))/F391)/(s/F337);
end
end

function Data = GetData
% E/N calibration table from ICCdFilterAnalsysisCode.m. Column 1 is E/N (Td),
% column 6 the N2+(B)/N2(C) ratio it is fitted against.
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
