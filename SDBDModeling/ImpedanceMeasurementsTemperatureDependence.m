%% HiokiImpedanceSweep.m
%
% Parse Hioki IM3536 impedance sweep files from a directory, extract the
% series-equivalent resistance, capacitance and loss tangent at every
% frequency point, and report how each quantity varies across files.
%
% FEATURES
%   - Automatic temperature extraction from filenames. Recognises patterns
%     like "25Degrees...", "T30_...", "30C_...", "_25deg_...", etc. Files
%     whose names contain no temperature are assigned NaN and grouped
%     separately as "unknown".
%   - Replicate averaging: files at the same temperature are combined.
%     The mean and standard deviation of every quantity (C, tan d, Re, Rp)
%     are computed per temperature, giving measurement uncertainty for free.
%   - Band-limited: points above MaxFreq (default 1 MHz) are discarded
%     before any analysis, per the IM3536's practical accuracy limits.
%
% USAGE
%   1. Set SweepDir to your directory, or leave it '' for a file browser.
%   2. Run the script. It produces:
%       - SweepTable:   struct array, one element per FILE
%       - TempTable:    struct array, one element per unique TEMPERATURE,
%                       holding the mean, std, and N of every quantity
%       - Summary figures with error bars where replicates exist
%
% DEPENDENCIES
%   natsort or natsortfiles (optional, for natural file ordering)
%
% NOTE ON UNITS
%   The Hioki reports |Z| in ohms and phase in degrees. Derived quantities:
%     Re(Z)  = |Z| cos(phase)          series resistance, ohm
%     Im(Z)  = |Z| sin(phase)          reactance, ohm (negative for a cap)
%     C_s    = -1 / (w Im(Z))          series capacitance, F
%     tan(d) = -Re(Z) / Im(Z)          loss tangent of the whole network
%     Rp     = 1 / (w C tan(d))        equivalent parallel resistance, ohm

%% Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SweepDir = '';            % '' = file browser
MaxFreq  = 1e6;          % Hz, discard above this
MinFreq  = 10;           % Hz, discard below this

% Frequencies at which to report trends
ReportFreqs = [100, 300, 1000, 3000, 10000, 30000, 100000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read directory
if isempty(SweepDir)
    SweepDir = uigetdir('E:\MATLAB\MATLAB\Data\2024\ImpedanceMEasurements\ImpedanceMEasurements\PET','Select directory of Hioki CSV files');
    if SweepDir == 0, error('No directory selected.'); end
end

listing = dir(fullfile(SweepDir,'*.csv'));
if isempty(listing), error('No CSV files found in %s',SweepDir); end

fnames = {listing.name};
try fnames = natsort(fnames); catch, fnames = sort(fnames); end
Nfiles = numel(fnames);
fprintf('Found %d CSV files in %s\n',Nfiles,SweepDir);

%% Parse each file and extract temperature
SweepTable = struct([]);

for n = 1:Nfiles
    fullpath = fullfile(SweepDir,fnames{n});
    [freq,Zmag,Phase] = ParseHioki(fullpath);

    % Band limit
    keep = (freq >= MinFreq) & (freq <= MaxFreq) & freq ~= 60 & freq ~=50 & freq ~= 70;
    freq = freq(keep); Zmag = Zmag(keep); Phase = Phase(keep);

    w  = 2*pi*freq;
    Re = Zmag .* cosd(Phase);
    Im = Zmag .* sind(Phase);
    Cs = -1./(w.*Im);
    td = -Re./Im;
    Rp = 1./(w.*Cs.*td);

    SweepTable(n).FileName = fnames{n};
    SweepTable(n).Temp     = ExtractTemperature(fnames{n});
    SweepTable(n).Freq     = freq;
    SweepTable(n).Zmag     = Zmag;
    SweepTable(n).Phase    = Phase;
    SweepTable(n).Re       = Re;
    SweepTable(n).Im       = Im;
    SweepTable(n).Cs       = Cs;
    SweepTable(n).TanDelta = td;
    SweepTable(n).Rp       = Rp;

    rband = (freq >= 100) & (freq <= 100e3);
    SweepTable(n).Cmean      = mean(Cs(rband));
    SweepTable(n).TanDelta1k = td(nearestIdx(freq,1000));

    if isnan(SweepTable(n).Temp)
        tstr = 'unknown';
    else
        tstr = sprintf('%g C',SweepTable(n).Temp);
    end
    fprintf('  %3d. %-35s  T = %s  (%d pts)\n',n,fnames{n},tstr,numel(freq));
end

%% Group by temperature, compute mean and std
allTemps = [SweepTable.Temp];

% Unique temperatures, NaN grouped as one bucket at the end
knownMask = ~isnan(allTemps);
uTemps = unique(allTemps(knownMask));
uTemps = uTemps(:);
if any(~knownMask)
    uTemps = [uTemps; NaN];
end
Ntemps = numel(uTemps);

% Frequency grid from the first file (all files should share it)
fgrid  = SweepTable(1).Freq;
Nfreqs = numel(fgrid);

TempTable = struct([]);

fprintf('\n=== Grouping by temperature ===\n');
for k = 1:Ntemps
    T = uTemps(k);
    if isnan(T)
        mask = isnan(allTemps);
        tlab = 'unknown';
    else
        mask = (allTemps == T);
        tlab = sprintf('%g C',T);
    end
    idx = find(mask);
    N = numel(idx);

    % Stack per-frequency vectors into matrices (Nreps x Nfreqs)
    tdMat = nan(N,Nfreqs);
    CsMat = nan(N,Nfreqs);
    ReMat = nan(N,Nfreqs);
    RpMat = nan(N,Nfreqs);
    CmeanVec = nan(N,1);

    for j = 1:N
        Nf = numel(SweepTable(idx(j)).Freq);
        Nuse = min(Nf,Nfreqs);
        tdMat(j,1:Nuse) = SweepTable(idx(j)).TanDelta(1:Nuse);
        CsMat(j,1:Nuse) = SweepTable(idx(j)).Cs(1:Nuse);
        ReMat(j,1:Nuse) = SweepTable(idx(j)).Re(1:Nuse);
        RpMat(j,1:Nuse) = SweepTable(idx(j)).Rp(1:Nuse);
        CmeanVec(j)     = SweepTable(idx(j)).Cmean;
    end

    TempTable(k).Temp          = T;
    TempTable(k).Label         = tlab;
    TempTable(k).Nreps         = N;
    TempTable(k).FileIdx       = idx;
    TempTable(k).Freq          = fgrid;
    TempTable(k).TanDelta_mean = mean(tdMat,1,'omitnan');
    TempTable(k).TanDelta_std  = std(tdMat,0,1,'omitnan');
    TempTable(k).Cs_mean       = mean(CsMat,1,'omitnan');
    TempTable(k).Cs_std        = std(CsMat,0,1,'omitnan');
    TempTable(k).Re_mean       = mean(ReMat,1,'omitnan');
    TempTable(k).Re_std        = std(ReMat,0,1,'omitnan');
    TempTable(k).Rp_mean       = mean(RpMat,1,'omitnan');
    TempTable(k).Rp_std        = std(RpMat,0,1,'omitnan');
    TempTable(k).Cmean_mean    = mean(CmeanVec,'omitnan');
    TempTable(k).Cmean_std     = std(CmeanVec,0,'omitnan');

    j1k = nearestIdx(fgrid,1000);
    if N > 1
        fprintf('  %8s: %d files  C = %.2f +- %.2f pF  td@1k = %.5f +- %.5f\n', ...
                tlab,N, ...
                TempTable(k).Cmean_mean*1e12,TempTable(k).Cmean_std*1e12, ...
                TempTable(k).TanDelta_mean(j1k),TempTable(k).TanDelta_std(j1k));
    else
        fprintf('  %8s: %d file   C = %.2f pF  td@1k = %.5f\n', ...
                tlab,N, ...
                TempTable(k).Cmean_mean*1e12, ...
                TempTable(k).TanDelta_mean(j1k));
    end
end

%% Print summary table with uncertainties
fprintf('\n%-8s %3s %14s','Temp','N','C (pF)');
for m = 1:numel(ReportFreqs)
    fprintf('  td@%-7s',FreqLabel(ReportFreqs(m)));
end
fprintf('\n');
for k = 1:Ntemps
    if TempTable(k).Nreps > 1
        cstr = sprintf('%.1f+-%.1f',TempTable(k).Cmean_mean*1e12,TempTable(k).Cmean_std*1e12);
    else
        cstr = sprintf('%.1f',TempTable(k).Cmean_mean*1e12);
    end
    fprintf('%-8s %3d %14s',TempTable(k).Label,TempTable(k).Nreps,cstr);
    for m = 1:numel(ReportFreqs)
        j = nearestIdx(fgrid,ReportFreqs(m));
        mu = TempTable(k).TanDelta_mean(j);
        sd = TempTable(k).TanDelta_std(j);
        if TempTable(k).Nreps > 1
            fprintf('  %.4f+-%.4f',mu,sd);
        else
            fprintf('  %.5f      ',mu);
        end
    end
    fprintf('\n');
end

%% Temperature coefficient fit (only for resolved temperatures)
resolvedIdx = find(~isnan([TempTable.Temp]));
if numel(resolvedIdx) >= 3
    Tvec = [TempTable(resolvedIdx).Temp];
    fprintf('\n=== Temperature coefficients ===\n');
    fprintf('%8s %14s %10s %12s %6s\n','f','d(td)/dT','td(25C)','TCC','R2');
    for m = 1:numel(ReportFreqs)
        j = nearestIdx(fgrid,ReportFreqs(m));
        tdvec = arrayfun(@(idx) TempTable(idx).TanDelta_mean(j),resolvedIdx);
        p = polyfit(Tvec,tdvec,1);
        pred = polyval(p,Tvec);
        ss_res = sum((tdvec-pred).^2);
        ss_tot = sum((tdvec-mean(tdvec)).^2);
        R2 = 1 - ss_res/max(ss_tot,eps);
        td25 = polyval(p,25);
        TCC = p(1)/td25*1e6;
        fprintf('%8s %+13.2fe-6 %10.5f %+11.0f ppm/K %6.3f\n', ...
                FreqLabel(ReportFreqs(m)),p(1)*1e6,td25,TCC,R2);
    end

    Cvec = [TempTable(resolvedIdx).Cmean_mean];
    p = polyfit(Tvec,Cvec,1);
    fprintf('\n  C(T) = %.2f %+.4f*T  pF\n',p(2)*1e12,p(1)*1e12);
    fprintf('  TCC(C) = %+.0f ppm/K\n',p(1)/mean(Cvec)*1e6);
end

%% Figures

% Colour map by temperature
if numel(resolvedIdx) >= 2
    cmap = parula(numel(resolvedIdx));
    cmapAll = zeros(Ntemps,3);
    for k = 1:numel(resolvedIdx)
        cmapAll(resolvedIdx(k),:) = cmap(k,:);
    end
    unknownIdx = find(isnan([TempTable.Temp]));
    for k = 1:numel(unknownIdx)
        cmapAll(unknownIdx(k),:) = [0.5 0.5 0.5];
    end
else
    cmapAll = lines(Ntemps);
end

% --- tan(d) vs frequency with error bands ---
figure('Name','TanDelta vs Frequency')
legs = {};
for k = 1:Ntemps
    f  = TempTable(k).Freq;
    mu = TempTable(k).TanDelta_mean;
    sd = TempTable(k).TanDelta_std;
    semilogx(f,mu,'Color',cmapAll(k,:),'LineWidth',1.6); hold on
    if TempTable(k).Nreps > 1
        fill([f(:); flipud(f(:))],[mu(:)+sd(:); flipud(mu(:)-sd(:))], ...
             cmapAll(k,:),'FaceAlpha',0.15,'EdgeColor','none');
    end
    legs{end+1} = sprintf('%s (N=%d)',TempTable(k).Label,TempTable(k).Nreps); %#ok<SAGROW>
end
xlabel('Frequency / Hz'); ylabel('tan \delta')
title('Loss Tangent vs Frequency'); grid on; set(gca,'XScale','log')
legend(legs,'Location','NorthWest'); set(gcf,'Visible','on')

% --- Capacitance vs frequency with error bands ---
figure('Name','Capacitance vs Frequency')
legs = {};
for k = 1:Ntemps
    f  = TempTable(k).Freq;
    mu = TempTable(k).Cs_mean*1e12;
    sd = TempTable(k).Cs_std*1e12;
    semilogx(f,mu,'Color',cmapAll(k,:),'LineWidth',1.6); hold on
    if TempTable(k).Nreps > 1
        fill([f(:); flipud(f(:))],[mu(:)+sd(:); flipud(mu(:)-sd(:))], ...
             cmapAll(k,:),'FaceAlpha',0.15,'EdgeColor','none');
    end
    legs{end+1} = TempTable(k).Label; %#ok<SAGROW>
end
xlabel('Frequency / Hz'); ylabel('Capacitance / pF')
title('Series Capacitance vs Frequency'); grid on; set(gca,'XScale','log')
legend(legs,'Location','NorthEast'); set(gcf,'Visible','on')

% --- tan(d) vs temperature at selected frequencies (with error bars) ---
if numel(resolvedIdx) >= 2
    figure('Name','TanDelta vs Temperature')
    plotFreqs = [100 1000 10000 100000];
    mkrs = {'o','s','d','^'};
    legs = {};
    for m = 1:numel(plotFreqs)
        j   = nearestIdx(fgrid,plotFreqs(m));
        Tv  = [TempTable(resolvedIdx).Temp];
        mu  = arrayfun(@(idx) TempTable(idx).TanDelta_mean(j),resolvedIdx);
        sd  = arrayfun(@(idx) TempTable(idx).TanDelta_std(j),resolvedIdx);
        errorbar(Tv,mu,sd,['-' mkrs{m}],'LineWidth',1.4,'MarkerSize',7, ...
                 'CapSize',8); hold on
        legs{end+1} = FreqLabel(plotFreqs(m)); %#ok<SAGROW>
    end
    xlabel('Temperature / °C'); ylabel('tan \delta')
    title('Loss Tangent vs Temperature')
    legend(legs,'Location','NorthEast'); grid on; set(gcf,'Visible','on')

    % --- Capacitance vs temperature (with error bars) ---
    figure('Name','Capacitance vs Temperature')
    Tv    = [TempTable(resolvedIdx).Temp];
    Cmean = [TempTable(resolvedIdx).Cmean_mean]*1e12;
    Cstd  = [TempTable(resolvedIdx).Cmean_std]*1e12;
    errorbar(Tv,Cmean,Cstd,'-ok','LineWidth',1.4,'MarkerSize',7,'CapSize',8);
    xlabel('Temperature / °C'); ylabel('Capacitance / pF')
    title('Capacitance vs Temperature'); grid on; set(gcf,'Visible','on')
end

% --- Parallel resistance vs frequency ---
figure('Name','Parallel Resistance vs Frequency')
sel = unique([1 max(round(Ntemps/2),1) Ntemps]);
legs = {};
for k = sel
    loglog(TempTable(k).Freq,TempTable(k).Rp_mean/1e6,'LineWidth',1.4); hold on
    legs{end+1} = TempTable(k).Label; %#ok<SAGROW>
end
xlabel('Frequency / Hz'); ylabel('R_p / M\Omega')
title('Equivalent Parallel Resistance')
legend(legs,'Location','NorthEast'); grid on; set(gcf,'Visible','on')

fprintf('\nDone. Results in SweepTable (per file) and TempTable (per temperature).\n');

%% ========================================================================
%  Local functions
%  ========================================================================

function [freq,Zmag,Phase] = ParseHioki(filepath)
% Parse a Hioki IM3536 CSV export file.
txt = fileread(filepath);
txt = strrep(txt,char(13),'');
ftok  = regexp(txt,'FREQ,([\d.E+-]+),Hz','tokens');
freq  = cellfun(@(c) str2double(c{1}),ftok(:));
ztok  = regexp(txt,'Z\[ohm\],PHASE\[deg\]\s*\n([\d.E+-]+),(-?[\d.]+)','tokens');
Zmag  = cellfun(@(c) str2double(c{1}),ztok(:));
Phase = cellfun(@(c) str2double(c{2}),ztok(:));
if numel(freq) ~= numel(Zmag)
    warning('ParseHioki: %s has %d freq / %d Z entries',filepath,numel(freq),numel(Zmag));
    N = min(numel(freq),numel(Zmag));
    freq = freq(1:N); Zmag = Zmag(1:N); Phase = Phase(1:N);
end
end

function T = ExtractTemperature(fname)
% Extract temperature in degrees C from a filename.
%
% Recognised patterns (case insensitive, tested in order):
%   25Degrees...          25_Degrees...      25.5Degrees...
%   T25_...               t25...             T_25...
%   25C_...               25c.csv            25C at end of stem
%   25deg...              25_deg...
%   Temp25...             temp_25...
%   _25_C_...
%
% Returns NaN if nothing matches or the extracted number is outside -40..300.
T = NaN;
patterns = {
    '(\d+\.?\d*)[\s_-]*[Dd]egrees?'          % 25Degrees, 25_degrees
    '[Tt][\s_-]*(\d+\.?\d*)'                 % T25, T_25
    '(\d+\.?\d*)[\s_-]*[Cc](?=[_.\s\d]|$)'   % 25C_, 25c.csv
    '(\d+\.?\d*)[\s_-]*[Dd]eg'               % 25deg
    '[Tt]emp[\s_-]*(\d+\.?\d*)'              % Temp25
};
for k = 1:numel(patterns)
    tok = regexp(fname,patterns{k},'tokens','once');
    if ~isempty(tok)
        val = str2double(tok{1});
        if val >= -40 && val <= 300
            T = val;
            return
        end
    end
end
end

function j = nearestIdx(vec,target)
[~,j] = min(abs(vec - target));
end

function s = FreqLabel(f)
if f >= 1e6,    s = sprintf('%.0fM',f/1e6);
elseif f >= 1e3, s = sprintf('%.0fk',f/1e3);
else,            s = sprintf('%.0f',f);
end
end