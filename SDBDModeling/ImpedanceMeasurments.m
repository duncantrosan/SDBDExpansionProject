%% HiokiImpedanceSweep.m
%
% Parse Hioki IM3536 impedance sweep files from a directory, extract the
% series-equivalent resistance, capacitance and loss tangent at every
% frequency point, and report how each quantity varies across files
% (temperature, time, etc.).
%
% USAGE
%   1. Set SweepDir to your directory of Hioki CSV exports, or leave it
%      empty to get a file browser.
%   2. Run the script. It produces:
%       - SweepTable: a struct array, one element per file, each holding
%         the full f/Z/phase/R/C/tand vectors
%       - Summary figures: C vs file, tan(d) vs file at selected
%         frequencies, and a full tan(d)-vs-f overlay
%
% The Hioki IM3536 format repeats a block of metadata for every frequency
% point. The parser pulls FREQ and Z/PHASE from each block using textscan,
% so it is insensitive to the number of header lines or the block layout
% changing between firmware versions.
%
% DEPENDENCIES
%   natsort or natsortfiles (optional, for natural file ordering)
%
% NOTE ON UNITS
%   The Hioki reports |Z| in ohms and phase in degrees. Everything below
%   is derived from those two numbers:
%     Re(Z)  = |Z| cos(phase)          series resistance, ohm
%     Im(Z)  = |Z| sin(phase)          reactance, ohm (negative for a cap)
%     C_s    = -1 / (w Im(Z))          series capacitance, F
%     tan(d) = -Re(Z) / Im(Z)          loss tangent of the WHOLE network
%
%   tan(d) = -Re/Im is the device's loss tangent directly. It equals the
%   ratio of real power dissipated to reactive power stored per cycle.
%   Physically, it is a resistance Rp in PARALLEL with the capacitor:
%
%           Rp = 1 / (w C tan(d))
%
%   At 1 kHz and 490 pF, tan(d) = 0.020 corresponds to Rp = 16.2 MOhm.
%   The series-equivalent resistance is Rs = tan(d)/(wC) = 6.5 kOhm.
%   Same physics, two representations; the parallel form is natural for a
%   leaky dielectric, the series form is what the sine-fit impedance code
%   returns.

%% Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SweepDir = '';            % set to '' for a file browser
MaxFreq  = 1e6;          % Hz, discard points above this (instrument limit)
MinFreq  = 10;           % Hz, discard points below this if desired

% Frequencies at which to report tan(d) trends across files
ReportFreqs = [100, 300, 1000, 3000, 10000, 30000, 100000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read directory
if isempty(SweepDir)
    SweepDir = uigetdir('E:\MATLAB\MATLAB\Data\2024\ImpedanceMEasurements\ImpedanceMEasurements\PET','Select directory of Hioki CSV files');
    if SweepDir == 0
        error('No directory selected.');
    end
end

listing = dir(fullfile(SweepDir,'*.csv'));
if isempty(listing)
    error('No CSV files found in %s',SweepDir);
end

% Natural sort if available, otherwise alphabetical (which is fine for
% timestamped Hioki filenames like 240709042801.csv)
fnames = {listing.name};
try
    fnames = natsort(fnames);
catch
    fnames = sort(fnames);
end
Nfiles = numel(fnames);
fprintf('Found %d CSV files in %s\n',Nfiles,SweepDir);

%% Parse each file
SweepTable = struct([]);

for n = 1:Nfiles
    fullpath = fullfile(SweepDir,fnames{n});
    [freq,Zmag,Phase] = ParseHioki(fullpath);

    % Band limit
    keep = (freq >= MinFreq) & (freq <= MaxFreq) & freq ~= 60 & freq ~=50 & freq ~= 70;
    freq  = freq(keep);
    Zmag  = Zmag(keep);
    Phase = Phase(keep);

    w  = 2*pi*freq;
    Re = Zmag .* cosd(Phase);       % Re(Z), series resistance
    Im = Zmag .* sind(Phase);       % Im(Z), negative for a capacitor
    Cs = -1./(w.*Im);               % series capacitance
    td = -Re./Im;                   % loss tangent

    % Parallel equivalent resistance: Rp = 1/(w*C*tan(d))
    Rp = 1./(w.*Cs.*td);

    SweepTable(n).FileName  = fnames{n};
    SweepTable(n).Freq      = freq;
    SweepTable(n).Zmag      = Zmag;
    SweepTable(n).Phase     = Phase;
    SweepTable(n).Re        = Re;
    SweepTable(n).Im        = Im;
    SweepTable(n).Cs        = Cs;
    SweepTable(n).TanDelta  = td;
    SweepTable(n).Rp        = Rp;

    % Scalar summaries for the reporting band
    rband = (freq >= 100) & (freq <= 100e3);
    SweepTable(n).Cmean     = mean(Cs(rband));
    SweepTable(n).TanDelta1k= td(nearestIdx(freq,1000));
end

%% Print summary table
fprintf('\n%-18s %9s','File','C(pF)');
for k = 1:numel(ReportFreqs)
    fprintf(' td@%6s',FreqLabel(ReportFreqs(k)));
end
fprintf('\n');
for n = 1:Nfiles
    fprintf('%-18s %9.2f',SweepTable(n).FileName,SweepTable(n).Cmean*1e12);
    for k = 1:numel(ReportFreqs)
        j = nearestIdx(SweepTable(n).Freq,ReportFreqs(k));
        fprintf(' %8.5f',SweepTable(n).TanDelta(j));
    end
    fprintf('\n');
end

% Change from first to last
fprintf('\nFirst -> Last file:\n');
for k = 1:numel(ReportFreqs)
    j1 = nearestIdx(SweepTable(1).Freq,ReportFreqs(k));
    jN = nearestIdx(SweepTable(Nfiles).Freq,ReportFreqs(k));
    a = SweepTable(1).TanDelta(j1);
    b = SweepTable(Nfiles).TanDelta(jN);
    fprintf('  tan(d) @%6s: %.5f -> %.5f  (%+.1f%%)\n', ...
            FreqLabel(ReportFreqs(k)),a,b,100*(b/a-1));
end
fprintf('  C: %.2f -> %.2f pF  (%+.2f%%)\n', ...
        SweepTable(1).Cmean*1e12,SweepTable(Nfiles).Cmean*1e12, ...
        100*(SweepTable(Nfiles).Cmean/SweepTable(1).Cmean-1));

%% Figures

% --- tan(d) vs frequency, all files overlaid ---
figure('Name','TanDelta vs Frequency')
cmap = parula(Nfiles);
for n = 1:Nfiles
    semilogx(SweepTable(n).Freq,SweepTable(n).TanDelta, ...
             'Color',cmap(n,:),'LineWidth',1.2); hold on
end
xlabel('Frequency / Hz')
ylabel('tan \delta')
title('Loss Tangent vs Frequency')
cb = colorbar; cb.Label.String = 'File index (time / temperature)';
clim([1 Nfiles])
grid on; set(gca,'XScale','log')
set(gcf,'Visible','on')

% --- Capacitance vs frequency, all files overlaid ---
figure('Name','Capacitance vs Frequency')
for n = 1:Nfiles
    semilogx(SweepTable(n).Freq,SweepTable(n).Cs*1e12, ...
             'Color',cmap(n,:),'LineWidth',1.2); hold on
end
xlabel('Frequency / Hz')
ylabel('Capacitance / pF')
title('Series Capacitance vs Frequency')
cb = colorbar; cb.Label.String = 'File index';
clim([1 Nfiles])
grid on; set(gca,'XScale','log')
set(gcf,'Visible','on')

% --- tan(d) at selected frequencies across files ---
figure('Name','TanDelta Trend Across Files')
plotFreqs = [100 1000 10000 100000];
markers = {'o','s','d','^'};
for k = 1:numel(plotFreqs)
    vals = nan(1,Nfiles);
    for n = 1:Nfiles
        j = nearestIdx(SweepTable(n).Freq,plotFreqs(k));
        vals(n) = SweepTable(n).TanDelta(j);
    end
    plot(1:Nfiles,vals,['-' markers{k}],'LineWidth',1.4,'MarkerSize',6); hold on
end
xlabel('File index (time / temperature order)')
ylabel('tan \delta')
title('Loss Tangent Trend Across Measurements')
legend(cellfun(@(f)sprintf('%s',FreqLabel(f)),num2cell(plotFreqs),'Uni',0), ...
       'Location','NorthEast')
grid on
set(gcf,'Visible','on')

% --- Parallel resistance vs frequency ---
figure('Name','Parallel Resistance vs Frequency')
for n = [1 round(Nfiles/2) Nfiles]
    loglog(SweepTable(n).Freq,SweepTable(n).Rp/1e6,'LineWidth',1.4); hold on
end
xlabel('Frequency / Hz')
ylabel('R_p / M\Omega')
title('Equivalent Parallel Resistance')
legend({SweepTable(1).FileName,SweepTable(round(Nfiles/2)).FileName, ...
        SweepTable(Nfiles).FileName},'Location','NorthEast')
grid on
set(gcf,'Visible','on')

fprintf('\nDone. Results are in the SweepTable struct array.\n');

%% ========================================================================
%  Local functions
%  ========================================================================

function [freq,Zmag,Phase] = ParseHioki(filepath)
% Parse a Hioki IM3536 CSV export file.
% Returns column vectors of frequency (Hz), |Z| (ohm), and phase (deg).
%
% The format repeats a metadata block for each frequency point:
%   No.NNN
%   FREQ,<value>,Hz
%   ...settings...
%   Z[ohm],PHASE[deg]
%   <Zmag>,<Phase>
%
% This parser reads the whole file as text and pulls the numbers out with
% regexp, which makes it insensitive to the exact number of settings lines.

txt = fileread(filepath);
% Handle both \r\n and \n line endings
txt = strrep(txt,char(13),'');

% Frequencies
ftok = regexp(txt,'FREQ,([\d.E+-]+),Hz','tokens');
freq = cellfun(@(c) str2double(c{1}),ftok(:));

% Z magnitude and phase, on the line after "Z[ohm],PHASE[deg]"
ztok = regexp(txt,'Z\[ohm\],PHASE\[deg\]\s*\n([\d.E+-]+),(-?[\d.]+)','tokens');
Zmag  = cellfun(@(c) str2double(c{1}),ztok(:));
Phase = cellfun(@(c) str2double(c{2}),ztok(:));

% Sanity check
if numel(freq) ~= numel(Zmag)
    warning('ParseHioki: %s has %d freq entries but %d Z/phase pairs', ...
            filepath,numel(freq),numel(Zmag));
    N = min(numel(freq),numel(Zmag));
    freq = freq(1:N); Zmag = Zmag(1:N); Phase = Phase(1:N);
end
end

function j = nearestIdx(vec,target)
% Index of the element in vec closest to target.
[~,j] = min(abs(vec - target));
end

function s = FreqLabel(f)
% Human-readable frequency label.
if f >= 1e6
    s = sprintf('%.0fM',f/1e6);
elseif f >= 1e3
    s = sprintf('%.0fk',f/1e3);
else
    s = sprintf('%.0f',f);
end
end