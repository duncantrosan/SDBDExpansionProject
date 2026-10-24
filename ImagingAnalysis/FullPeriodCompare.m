%% Delayed vs prompt sides over a full period
%
% The echelle spectrum (EchelleFilterFactors.m) is integrated over whole
% periods, so its N2+ share of the 391 filter signal (FNSFraction, ~0.45) is
% only strictly valid for light integrated the same way. This script redoes
% the per-side N2(C)/N2+(B) vs ignition-delay comparison on full-period light,
% from two
% independent sources:
%   Sweep  - the cached 2 us gate sweep (BigHexDerived.mat), summed over every
%            gate of the period, per side
%   Image  - the separate full-period exposures in BigHex/FullPeriodImaging
% plus the positive half-wave alone (summed over its gates) as the reference.
%
% Two corrections for the N2 SPS light leaking into the 391 filter:
%   Fraction - N2+ = FNSFraction * S391 for every side (the established method)
%   Subtract - N2+ = S391 - Kappa * S337, with Kappa set so that the
%              electrode-integrated N2+ share over the full period equals the
%              echelle's FNSFraction. Same electrode-wide value as 'Fraction',
%              but a side's lower raw ratio is read as less N2+, which is what
%              the leakage physics says (the leak scales with S337, not S391).
%
% Ignition delays (which sides are late) come from the time-resolved sweep
% (SideTable in the cache), per half-wave.
%
% Inputs: Output/Cache (run ICCDFilterRatioEoN.m first) and the full-period CSVs.
% Outputs: Output/FullPeriod. Figures are left open and saved as .png and .fig.
%
% ASSUMPTIONS
%   - The echelle spectrum is from the Kapton small-hex electrode, not Big Hex;
%     its N2+ share is taken to hold for Big Hex.
%   - The 391 full-period exposure has no background of its own: the 391 sweep
%     background is used, and any remaining offset is removed using pixels well
%     away from the electrode.

clear

%% Section 1 - Options
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
CacheDir = fullfile(RepoRoot,'Output','Cache');
OutDir   = fullfile(RepoRoot,'Output','FullPeriod');
BigHex   = fullfile(RepoRoot,'Data','FullImagingDataSet','Fulldata','12_19_23Meanenergy','BigHex');
FP337File = fullfile(BigHex,'FullPeriodImaging','FullPeriodData','337Line','BigHexFullPeriod337Line 2023.csv');
FP391File = fullfile(BigHex,'FullPeriodImaging','FullPeriodData','391Line','BigHexFullPeriod391Line 2023 December 19 17_15_39.csv');
BG337Glob = fullfile(BigHex,'FullPeriodImaging','BigHexFullPeriod337LineBG*.csv');
BG391File = fullfile(BigHex,'SweepOverTime','3912usSweep','BigHex2ussweep391LineBG 2023 December 19 17_13_38.csv');
SpecFile  = fullfile(RepoRoot,'Data','ICCD_ImagingData','Echelledata','AbsoultelyCalibratedtrial', ...
                     'Absoultecalibrationtest','Test','KaptonSmallHex.spa');
TransFile = fullfile(RepoRoot,'Data','Stuff','FilterTransmitance.csv');
TreatFile = fullfile(RepoRoot,'BiggHexImaging','ImageTreatment.mat');

DelayThreshold = 4;                    % us after the median side, as before
DelayBinEdges  = [-Inf 0 4 8 12 16 24 Inf];
DelayBinMid    = [-4 2 6 10 14 20 30];
DarkFrac       = 0.01;                 % "dark" = summed sweep 337 below this share of its max
                                       % (scattered light keeps the whole frame at 0.2-1%)
SaturationCounts = 60000;
ReportEoN      = 0;                    % 1 = also report E/N (Td) next to the ratio

if ~exist(OutDir,'dir'), mkdir(OutDir); end

%% Section 2 - Cached sweep data and filter factors
D  = load(fullfile(CacheDir,'BigHexDerived.mat'));
Gm = load(fullfile(CacheDir,'HexagonGeometryRefined.mat'),'RectLabel','SideMid','SideCell','SideNum');
if ~strcmpi(D.ContaminationMode,'fraction')
    error('The cache was built with ContaminationMode = ''%s''; rebuild it with ''fraction''.',D.ContaminationMode);
end
FF = EchelleFilterFactors(SpecFile,TransFile,0,1);
f    = FF.FNSFraction;
T337 = D.FilterFactor337;
T391 = D.FilterFactor391;
EchelleRaw = (FF.SPS391 + FF.FNS391)/FF.S337;
Nsides = max(Gm.RectLabel(:));

T = D.SideTable;
DelayPos = T.Delay_us(strcmp(T.HalfWave,'Positive'));
DelayNeg = T.Delay_us(strcmp(T.HalfWave,'Negative'));

% Per-side raw sums: 337, and the measured 391 (the cache holds the anchored
% N2+ signal, AnchorGain * FNSFraction * S391)
Raw391 = D.SideFNS/(f*D.AnchorGain);
gAll = [];  for h = 1:numel(D.HumpStart), gAll = [gAll D.HumpStart(h):D.HumpEnd(h)]; end %#ok<AGROW>
gPos = []; gNeg = [];
for h = 1:numel(D.HumpStart)
    g = D.HumpStart(h):D.HumpEnd(h);
    if D.HalfWave(g(1)) > 0, gPos = g; else, gNeg = g; end
end
Src.Name = {'Positive half-wave (sweep)','Full period (sweep summed)','Full period (exposure)'};
Src.A = {sum(D.Side337(:,gPos),2), sum(D.Side337(:,gAll),2), []};
Src.B = {sum(Raw391(:,gPos),2),    sum(Raw391(:,gAll),2),    []};

%% Section 3 - Full-period exposures, on the same pixels as the sweep
IT = load(TreatFile);
CropRect = IT.ImageHandeling.Crop;
Stk = load(fullfile(CacheDir,'BigHexFilterStacks.mat'),'S337');
Ref = double(sum(Stk.S337,3));
clear Stk

BGf = dir(BG337Glob);
BGs = zeros(1024,1024,numel(BGf));
for k = 1:numel(BGf), BGs(:,:,k) = ReadPimaxCSV(fullfile(BGf(k).folder,BGf(k).name)); end
BG337 = median(BGs,3);                      % robust to the odd frame with stray light
BG391 = ReadPimaxCSV(BG391File);
Raw337 = ReadPimaxCSV(FP337File);
Raw391i = ReadPimaxCSV(FP391File);
fprintf('\nFull-period exposures: max raw counts 337 = %.0f, 391 = %.0f\n',max(Raw337(:)),max(Raw391i(:)));
if max(Raw337(:)) >= SaturationCounts || max(Raw391i(:)) >= SaturationCounts
    warning('Full-period exposure may be saturated (>= %d counts).',SaturationCounts);
end
I337 = medfilt2(imcrop(Raw337 - BG337,CropRect),[3 3]);
I391 = medfilt2(imcrop(Raw391i - BG391,CropRect),[3 3]);

% Align both exposures to the sweep (the camera may have moved slightly)
R = imref2d(size(Ref));
t337 = imregcorr(I337/max(I337(:)),Ref/max(Ref(:)),'translation');
t391 = imregcorr(I391/max(I391(:)),Ref/max(Ref(:)),'translation');
Inside = imwarp(true(size(I337)),t337,'OutputView',R,'FillValues',0) & ...
         imwarp(true(size(I391)),t391,'OutputView',R,'FillValues',0);
I337 = imwarp(I337,t337,'OutputView',R,'FillValues',0);
I391 = imwarp(I391,t391,'OutputView',R,'FillValues',0);
fprintf('Registration to the sweep: 337 shift (%.1f, %.1f) px, 391 shift (%.1f, %.1f) px\n', ...
        t337.T(3,1:2),t391.T(3,1:2));

% Remaining offset from pixels far from any plasma
Dark = imerode(Ref < DarkFrac*max(Ref(:)),strel('disk',15)) & Inside;
Off337 = median(I337(Dark));  Off391 = median(I391(Dark));
I337 = I337 - Off337;  I391 = I391 - Off391;
fprintf('Residual dark offsets removed: 337 %.1f, 391 %.1f counts/px (%d dark px)\n',Off337,Off391,nnz(Dark));

rm = Gm.RectLabel > 0;
Src.A{3} = accumarray(Gm.RectLabel(rm),I337(rm),[Nsides 1]);
Src.B{3} = accumarray(Gm.RectLabel(rm),I391(rm),[Nsides 1]);

%% Section 4 - Ratios and the delayed vs prompt comparison
% Anchoring: the sweep and the full-period exposures were not taken with the
% same 337/391 settings (accumulations, gain), and the ICCD and echelle
% responses differ. Over a full period both instruments see the same light, so
% each source is scaled so that its whole-frame, full-period 391/337 equals the
% echelle's. The sweep uses the factor from ICCDFilterRatioEoN.m (cached), the
% exposures their own. After anchoring, 'Fraction' uses the echelle N2+ share
% and 'Subtract' the echelle Kappa directly.
Gain = [D.AnchorGain D.AnchorGain EchelleRaw/(sum(I391(Inside))/sum(I337(Inside)))];
fprintf('\nEchelle full-period raw 391/337 = %.3f; anchor x%.2f (sweep), x%.2f (exposures)\n', ...
        EchelleRaw,Gain(1),Gain(3));
Corr = {'Fraction','Subtract'};
Res  = table;
CBSide = cell(3,2);
for s = 1:3
    A = Src.A{s};
    B = Gain(s)*Src.B{s};
    for c = 1:2
        if c == 1, N2p = f*B; else, N2p = B - FF.Kappa*A; end
        cb = (A/T337)./(N2p/T391);              % N2(C)/N2+(B) photons
        cb(N2p <= 0) = NaN;                      % no N2+ left after the leakage
        CBSide{s,c} = cb;
        for p = 1:2
            if p == 1, Dl = DelayPos; PolName = 'Positive'; else, Dl = DelayNeg; PolName = 'Negative'; end
            ok = ~isnan(Dl) & ~isnan(cb);
            dl = ok & Dl >= DelayThreshold;  pr = ok & Dl < DelayThreshold;
            if nnz(dl) < 3 || nnz(pr) < 3, continue, end
            [rho,prho] = corr(Dl(ok),cb(ok),'type','Spearman');
            Row = table(Src.Name(s),Corr(c),{PolName},nnz(N2p <= 0),median(cb(pr)),median(cb(dl)), ...
                  100*(median(cb(dl))/median(cb(pr)) - 1),ranksum(cb(dl),cb(pr)),rho,prho, ...
                  'VariableNames',{'Source','Correction','DelayFrom','NoN2plus','CBPrompt','CBDelayed', ...
                  'Diff_pct','RanksumP','SpearmanDelay','SpearmanP'});
            if ReportEoN == 1
                Row.EoNPrompt_Td  = polyval(D.EoNCoeff,1/Row.CBPrompt);
                Row.EoNDelayed_Td = polyval(D.EoNCoeff,1/Row.CBDelayed);
            end
            Res = [Res; Row]; %#ok<AGROW>
        end
    end
end
disp(Res(strcmp(Res.DelayFrom,'Positive'),:))
writetable(Res,fullfile(OutDir,'FullPeriodDelayedVsPrompt.csv'));
writetable(table((1:Nsides)',Gm.SideCell,Gm.SideNum,DelayPos,DelayNeg, ...
    CBSide{1,1},CBSide{2,1},CBSide{3,1},CBSide{1,2},CBSide{2,2},CBSide{3,2}, ...
    'VariableNames',{'Side','Cell','SideInCell','DelayPos_us','DelayNeg_us', ...
    'CBPosHalf_Fraction','CBSweepFull_Fraction','CBExposure_Fraction', ...
    'CBPosHalf_Subtract','CBSweepFull_Subtract','CBExposure_Subtract'}), ...
    fullfile(OutDir,'FullPeriodSideRatio_EchelleAnchored.csv'));

%% Section 5 - Figures
YLab = 'N_2(C) / N_2^+(B)';
Col = [0.55 0.55 0.55; 0.15 0.35 0.85; 0.85 0.25 0.1];
fig = figure('Color','w','Position',[40 40 1500 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for p = 1:2
    if p == 1, Dl = DelayPos; PolName = 'positive'; else, Dl = DelayNeg; PolName = 'negative'; end
    for c = 1:2
        nexttile; hold on
        Bin = discretize(Dl,DelayBinEdges);
        for s = 1:3
            V = CBSide{s,c};
            ok = ~isnan(Bin) & ~isnan(V);
            m  = accumarray(Bin(ok),V(ok),[numel(DelayBinMid) 1],@mean,NaN);
            se = accumarray(Bin(ok),V(ok),[numel(DelayBinMid) 1],@(x) std(x)/sqrt(numel(x)),NaN);
            errorbar(DelayBinMid,m,se,'o-','Color',Col(s,:),'MarkerFaceColor',Col(s,:),'LineWidth',1.3)
        end
        xline(DelayThreshold,'k--','HandleVisibility','off'); box on
        xlabel(sprintf('Side ignition delay in the %s half-wave / \\mus',PolName))
        ylabel([YLab ' (mean \pm s.e.m.)'])
        title({sprintf('%s correction, echelle-anchored',Corr{c}),sprintf('delay from the %s half-wave',PolName)})
        if p == 1 && c == 1, legend(Src.Name,'Location','northwest'); end
    end
    nexttile
    if p == 1
        % Do the two full-period sources agree side by side?
        scatter(CBSide{2,1},CBSide{3,1},18,DelayPos,'filled'); hold on
        lims = [min([CBSide{2,1};CBSide{3,1}]) max([CBSide{2,1};CBSide{3,1}])];
        plot(lims,lims,'k:'); axis equal tight; box on
        c = colorbar; c.Label.String = 'Positive-half-wave delay / \mus';
        xlabel(['Full period, sweep summed: ' YLab]); ylabel(['Full period, exposure: ' YLab])
        rr = corr(CBSide{2,1},CBSide{3,1},'rows','complete');
        title({'Per-side full-period ratio (fraction, anchored)',sprintf('sweep vs exposure, r = %.2f',rr)})
    else
        M = nan(size(Gm.RectLabel));
        M(rm) = CBSide{1,1}(Gm.RectLabel(rm));
        imagesc(M,'AlphaData',~isnan(M)); axis image off; colormap(gca,flipud(parula))
        c = colorbar; c.Label.String = ['Positive half-wave ' YLab];
        hold on
        k = find(DelayPos >= DelayThreshold);
        plot(Gm.SideMid(k,1),Gm.SideMid(k,2),'mo','MarkerSize',6,'LineWidth',1.3)
        title({'Positive half-wave ratio per side, anchored','circles: late in the positive half-wave'})
    end
end
SaveFigure(fig,OutDir,'FullPeriodDelayedVsPrompt')

%% Local functions
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

function SaveFigure(fig,Folder,Name)
% Save a figure as .png and .fig and leave it open, titled with its file name.
fig.Name = Name;
fig.Visible = 'on';
exportgraphics(fig,fullfile(Folder,[Name '.png']))
savefig(fig,fullfile(Folder,[Name '.fig']))
end
