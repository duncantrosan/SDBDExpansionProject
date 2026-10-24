function F = EchelleFilterFactors(SpecFile,TransFile,MakePlot,PhotonUnits)
%% EchelleFilterFactors
% Derives the 337/391 bandpass-filter correction factors from a measured
% echelle spectrum instead of the hard coded peak transmissions
% (FilterFactor337 = 0.2815, FilterFactor391 = 0.4134) and the fixed
% N2+ fraction of the 391 filter signal (Factor391Percentage = 0.4533).
%
% The 391 filter (FWHM ~10 nm) passes the N2+ FNS (0-0) 391.4 nm and
% (1-1) 388.4 nm bands, but also the N2 SPS Delta v = -3 sequence
% ((4-7) 385.8, (3-6) 389.5, (2-5) 394.3, (1-4) 399.8, (0-3) 405.9 nm).
% That SPS leakage comes from the same N2(C) state as the 337 nm band, so it
% scales with the 337 signal rather than being a fixed share of the 391
% signal:
%
%     S391_FNS = S391 - Kappa * S337
%
% Kappa (SPS counts through the 391 filter per count through the 337 filter)
% depends on the N2(C) vibrational distribution but not on E/N, so it can be
% taken from one time-integrated spectrum and applied to time-resolved
% images. A fixed fraction cannot: the FNS share itself changes with E/N.
%
% INPUTS
%   SpecFile    - echelle .spa export (header line, then wavelength/intensity)
%   TransFile   - FilterTransmitance.csv (col 1 wavelength, col 4 337 filter,
%                 col 8 390 filter; fractions, not percent)
%   MakePlot    - optional, 1 to plot the spectrum with the band windows
%   PhotonUnits - optional, 1 (default) if the .spa is in energy units and
%                 should be converted to photon flux (x lambda) to match an
%                 ICCD; 0 if it is already in counts
%
% OUTPUT struct F
%   T337eff, T391eff - transmission averaged over the target band
%   Kappa            - SPS leakage coefficient (counts/counts)
%   FNSFraction      - N2+ share of the 391 filter signal for THIS spectrum
%   Band             - background-subtracted band integrals, for checking
%
% METHOD NOTES
%   - The echelle noise floor is a population of positive single-pixel
%     spikes, mean 15-50 counts/point, which is a large share of the weak
%     391-region bands. It is removed by subtracting the MEAN of flanking
%     band-free windows (mean, not median, so the spikes inside each band
%     window cancel on average).
%   - SPS (3-6) sits under the FNS band and cannot be measured directly. It is
%     taken as the geometric mean of the clean (2-5) and (4-7) bands from the
%     same Delta v = -3 sequence.
%   - ICCD photocathode QE is taken as equal at 337 and 391 nm.

if nargin < 3, MakePlot = 0; end
if nargin < 4, PhotonUnits = 1; end

%% Read the spectrum
fid = fopen(SpecFile,'r');
C = textscan(fid,'%f%f','HeaderLines',1);
fclose(fid);
wl = C{1};
Iraw = C{2};
I = Iraw;
if PhotonUnits == 1
    I = I.*wl/337;    % relative photon flux, normalised at 337 nm
end

%% Read the filter transmission (skips the multi-line comment in the file)
Lines = readlines(TransFile);
T = nan(numel(Lines),3);
for k = 2:numel(Lines)
    v = str2double(split(Lines(k),','));
    if numel(v) >= 8 && all(isfinite(v([1 4 8])))
        T(k,:) = v([1 4 8]);
    end
end
T = T(all(isfinite(T),2),:);
T(:,2:3) = max(T(:,2:3),0);
T337 = interp1(T(:,1),T(:,2),wl,'linear',0);
T391 = interp1(T(:,1),T(:,3),wl,'linear',0);

%% Windows (nm). SPS and FNS bands are degraded to the violet, so each window
% runs from below the band up to just past its head.
W.SPS337 = [329.5 337.25];   % (0-0) 337.1, (1-1) 333.9, (2-2) 331.0
W.Pass337= [326.0 348.0];    % where T337 > ~1%
W.FNS    = [385.9 391.50];   % N2+ (0-0) 391.4 and (1-1) 388.4, plus SPS (3-6)
W.SPS47  = [383.5 385.85];
W.SPS25  = [391.6 394.35];
W.SPS14  = [396.5 399.90];
W.SPS03  = [403.0 406.00];
Bg337Win = {[322 328],[338 343]};
Bg391Win = {[382.0 383.3],[395.0 396.3]};

In    = @(w) wl >= w(1) & wl <= w(2);
Integ = @(y,w) trapz(wl(In(w)),y(In(w)));
BgMean = @(ws) mean(I(In(ws{1}) | In(ws{2})));

Bg337 = BgMean(Bg337Win);
Bg391 = BgMean(Bg391Win);
I337 = I - Bg337;
I391 = I - Bg391;

%% Effective transmissions over the target bands
F.T337eff = Integ(T337.*I337,W.SPS337)/Integ(I337,W.SPS337);
F.T391eff = Integ(T391.*I391,W.FNS)/Integ(I391,W.FNS);

%% Band integrals (no filter) and the SPS (3-6) estimate
B.SPS47 = Integ(I391,W.SPS47);
B.SPS25 = Integ(I391,W.SPS25);
B.SPS14 = Integ(I391,W.SPS14);
B.SPS03 = Integ(I391,W.SPS03);
B.SPS36 = sqrt(max(B.SPS25,0)*max(B.SPS47,0));
B.FNSwin= Integ(I391,W.FNS);
B.FNS   = B.FNSwin - B.SPS36;

% (3-6) has the (2-5) shape shifted 4.84 nm to the blue; weight it by the
% filter transmission over that shifted window
T391at36 = mean(T391(In(W.SPS25 - 4.84)));

%% Signals through each filter
SPS391 = Integ(T391.*I391,W.SPS47) + Integ(T391.*I391,W.SPS25) + ...
         Integ(T391.*I391,W.SPS14) + Integ(T391.*I391,W.SPS03) + ...
         B.SPS36*T391at36;
FNS391 = Integ(T391.*I391,W.FNS) - B.SPS36*T391at36;
S337   = Integ(T337.*I337,W.Pass337);

F.Kappa       = SPS391/S337;
F.FNSFraction = FNS391/(FNS391 + SPS391);
F.SPS391      = SPS391;
F.FNS391      = FNS391;
F.S337        = S337;
F.Bg337       = Bg337;
F.Bg391       = Bg391;
F.Band        = B;
F.Windows     = W;
F.PhotonUnits = PhotonUnits;

fprintf('\nEchelle filter factors (%s)\n',SpecFile);
fprintf('  Background, 337 / 391 region : %.1f / %.1f counts per point\n',Bg337,Bg391);
fprintf('  SPS Dv=-3 bands (0-3 1-4 2-5 3-6* 4-7): %.0f %.0f %.0f %.0f %.0f  (*estimated)\n', ...
        B.SPS03,B.SPS14,B.SPS25,B.SPS36,B.SPS47);
fprintf('  T337 effective            : %.4f   (old FilterFactor337 = 0.2815)\n',F.T337eff);
fprintf('  T391 effective            : %.4f   (old FilterFactor391 = 0.4134)\n',F.T391eff);
fprintf('  Kappa (SPS leak 391/337)  : %.4f\n',F.Kappa);
fprintf('  FNS fraction of 391 signal: %.3f   (old Factor391Percentage = 0.4533)\n',F.FNSFraction);

if MakePlot == 1
    figure
    tiledlayout(2,1)
    nexttile
    plot(wl,I337,'k'); hold on
    ylim([-200 max(I337(In(W.SPS337)))*1.05]); ylabel('Intensity - background')
    yyaxis right
    plot(wl,T337,'b','LineWidth',1.5); ylabel('Transmission')
    xlim([320 362]); xlabel('Wavelength / nm'); title('337 filter')
    nexttile
    plot(wl,I391,'k'); hold on
    yl = [-200 max(I391(In([380 407])))*1.05];
    ylim(yl); ylabel('Intensity - background')
    Names = {'SPS47','FNS','SPS25','SPS14','SPS03'};
    for k = 1:numel(Names)
        w = W.(Names{k});
        patch([w(1) w(2) w(2) w(1)],[yl(1) yl(1) yl(2) yl(2)],'r', ...
              'FaceAlpha',0.06,'EdgeColor','none');
        text(mean(w),0.93*yl(2),Names{k},'HorizontalAlignment','center','FontSize',8)
    end
    yyaxis right
    plot(wl,T391,'r','LineWidth',1.5); ylabel('Transmission')
    xlim([378 408]); xlabel('Wavelength / nm'); title('391 filter')
end
end
