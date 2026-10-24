%% Fig 13 - old vs new
%
% Rebuilds paper Fig 13 (plasma charge derivative dQp/dUr against the
% electrode voltage, overlaid with four imaging quantities) twice, side by side:
%   Old - stamped hexagons + 1.7 us trigger offset (PerimeterCodev2.m as it was)
%   New - refitted hexagons + emission-aligned timing (ICCDFilterRatioEoN.m)
%
% Panels, as in the paper:
%   A  open perimeter ratio zeta
%   B  plasma area from the thresholded whole-electrode image
%   C  plasma area estimated from the per-side plasma extent with Eq 14
%   D  total 337 intensity inside the plasma mask
%
% Reads only BigHexDerived.mat (written by ICCDFilterRatioEoN.m, Section 10)
% and the two scope traces, so it runs in seconds. Rerun ICCDFilterRatioEoN.m
% only when the image processing settings change.
%
% Everything is plotted against Ur = V - I*Resr (the paper's U), for both old
% and new, so old vs new differs only in timing and hexagon geometry.

clear
close all

%% Section 1 - Options
RepoRoot  = fileparts(fileparts(mfilename('fullpath')));
OutDir    = fullfile(RepoRoot,'Data','FullImagingDataSet','Results','EoN_BigHex');
D = load(fullfile(OutDir,'BigHexDerived.mat'));
ScopeOn   = D.ScopeFile;
ScopeOff  = fullfile(fileparts(fileparts(ScopeOn)),'PlasmaOff','bigkex1vpp3000hz_ALL.csv');

PhaseBins        = 1000;   % bins per drive cycle
SGOrder          = 4;      % paper: fourth-degree Savitzky-Golay derivative
SGWindow         = 41;
TurningPointFrac = 0.10;   % paper: mask where dUr/dphase < 10% of its maximum

%% Section 2 - Off measurement: Resr and Coff
[tOff,VOff,IOff] = ReadScope(ScopeOff,D.AmpGain,D.CurrentGain/1000);
fOff = FitFrequency(tOff,VOff,D.Fdrive);
w    = 2*pi*fOff;
Z    = Phasor(tOff,VOff,fOff)/Phasor(tOff,IOff,fOff);
Resr = real(Z);
Coff = -1/(w*imag(Z));
fprintf('Off measurement at %.1f Hz: Resr = %.2f kOhm, Coff = %.1f pF\n',fOff,Resr/1e3,Coff*1e12);

%% Section 3 - On measurement: phase-averaged Ur, Qp and dQp/dUr
[t,V,I] = ReadScope(ScopeOn,D.AmpGain,D.CurrentGain/1000);
f  = FitFrequency(t,V,D.Fdrive);
P  = 1/f;
dt = median(diff(t));
Ncyc = floor((t(end) - t(1))/P);
k  = t < t(1) + Ncyc*P;                      % whole cycles only
t = t(k);  V = V(k);  I = I(k) - mean(I(k));
Ur = V - I*Resr;
Q  = cumtrapz(t,I);  Q = Q - mean(Q);
Qp = Q - Coff*Ur;                            % paper Eqs 7-8

Phase = mod((t - t(1))/P,1);
b   = min(floor(Phase*PhaseBins) + 1,PhaseBins);
cnt = accumarray(b,1,[PhaseBins 1]);
Ub  = accumarray(b,Ur,[PhaseBins 1])./cnt;
Qb  = accumarray(b,Qp,[PhaseBins 1])./cnt;
Vb  = accumarray(b,V, [PhaseBins 1])./cnt;
[~,g] = sgolay(SGOrder,SGWindow);
dU = CircDeriv(Ub,g);
dQ = CircDeriv(Qb,g);
dQdU = dQ./dU;
dQdU(abs(dU) < TurningPointFrac*max(abs(dU))) = NaN;
PhaseC = ((0:PhaseBins-1)' + 0.5)/PhaseBins;
fprintf('On measurement at %.1f Hz, %d cycles: Ur amplitude %.0f V\n',f,Ncyc,max(abs(Ub)));

% Ur at each gate. Phase is referenced to the first scope sample.
UrAt = @(tg) interp1([PhaseC-1; PhaseC; PhaseC+1],[Ub; Ub; Ub],mod((tg - t(1))/P,1));
tNew = (D.GateDelay + D.GateWidth/2)*1e-6 + D.ScopeTimeShift;
tOld = D.GateDelay*1e-6 - D.OldTriggerOffset;
UrNew = UrAt(tNew);
UrOld = UrAt(tOld);

%% Section 4 - Imaging quantities, old and new
% Eq 14: Ap = Ncells*(6*zp*s - 2*sqrt(3)*zp^2) + Pout*zp, with zp the mean
% plasma extent of the sides that are lit (the equation assumes no partial
% discharge) and s the hexagon side length in the same image scale.
Eq14 = @(Ext,s) D.Ncells*(6*Ext*s - 2*sqrt(3)*Ext.^2) + D.OuterArcMM*Ext;
ExtOld = MeanLitExtent(D.SideExtent(:,:,1));
ExtNew = MeanLitExtent(D.SideExtent(:,:,2));

Img.Name  = {'Open perimeter ratio \zeta','Plasma area (mask) / cm^2', ...
           'Estimated plasma area (Eq 14) / cm^2','Total intensity / 10^9 counts'};
Img.Old   = {D.ZetaOld, D.PlasmaArea/100, Eq14(ExtOld,D.SideLengthOldMM)/100, D.TotalInt/1e9};
Img.New   = {D.Zeta,    D.PlasmaArea/100, Eq14(ExtNew,D.SideLengthMM)/100,    D.TotalInt/1e9};

%% Section 5 - Plot: 4 rows (panels A-D) x 2 columns (old, new)
HWNames  = {'Negative','Positive'};
Neg      = dU < 0;                           % negative half-wave: Ur falling
Col      = [0.15 0.3 0.8; 0.85 0.3 0.1];
fig = figure('Color','w','Position',[50 30 1250 1300]);
tl = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
Version = {'Old: stamped hexagons, 1.7 \mus offset', ...
           sprintf('New: refitted hexagons, %.0f \\mus alignment',D.ScopeTimeShift*1e6)};
for p = 1:4
    for c = 1:2
        nexttile
        yyaxis left
        plot(Ub(Neg),dQdU(Neg)*1e12,'-','Color',[0.4 0.4 0.4],'LineWidth',1.2); hold on
        plot(Ub(~Neg),dQdU(~Neg)*1e12,'--','Color',[0.4 0.4 0.4],'LineWidth',1.2)
        ylabel('dQ_p/dU_r / pF'); ylim([-50 max(dQdU*1e12,[],'omitnan')*1.1])
        set(gca,'YColor',[0.3 0.3 0.3])
        yyaxis right
        if c == 1, Y = Img.Old{p}; U = UrOld; else, Y = Img.New{p}; U = UrNew; end
        for h = 1:numel(D.HumpStart)
            gg = D.HumpStart(h):D.HumpEnd(h);
            hw = 1 + (D.HalfWave(gg(1)) > 0);
            plot(U(gg),Y(gg),'o-','Color',Col(hw,:),'MarkerSize',3,'LineWidth',1.2)
        end
        ylabel(Img.Name{p}); set(gca,'YColor','k')
        xlabel('U_r / V'); xlim([-1.1 1.1]*max(abs(Ub))); box on
        if p == 1, title(Version{c}); end
        ax = gca;  ax.Toolbar.Visible = 'off';   % keep it out of the exported image
        text(0.02,0.92,char('A' + p - 1),'Units','normalized','FontWeight','bold','FontSize',14)
    end
end
lg = legend({'dQ_p/dU_r, negative','dQ_p/dU_r, positive','Imaging, negative','Imaging, positive'}, ...
            'Orientation','horizontal');
lg.Layout.Tile = 'south';
exportgraphics(fig,fullfile(OutDir,'Fig13_OldVsNew.png'),'Resolution',200)

%% Section 6 - Numbers and Origin export
for h = 1:numel(D.HumpStart)
    gg = D.HumpStart(h):D.HumpEnd(h);
    hw = 1 + (D.HalfWave(gg(1)) > 0);
    for c = 1:2
        if c == 1, Z_ = D.ZetaOld(gg); U = UrOld(gg); Name = 'old'; else, Z_ = D.Zeta(gg); U = UrNew(gg); Name = 'new'; end
        i0 = find(Z_ > 0.01,1);  i9 = find(Z_ >= 0.9,1);
        if isempty(i9)
            fprintf('%s half-wave, %s: first light at Ur = %.0f V, zeta never reaches 0.9 (max %.2f)\n', ...
                    HWNames{hw},Name,U(i0),max(Z_));
        else
            fprintf('%s half-wave, %s: first light at Ur = %.0f V, zeta >= 0.9 at %.0f V (%.0f V later), max %.2f\n', ...
                    HWNames{hw},Name,U(i0),U(i9),abs(U(i9)-U(i0)),max(Z_));
        end
    end
end
writetable(table(D.GateDelay(:),D.HalfWave(:),UrOld(:),UrNew(:),D.ZetaOld(:),D.Zeta(:), ...
                 D.PlasmaArea(:)/100,Img.Old{3}(:),Img.New{3}(:),ExtOld(:),ExtNew(:),D.TotalInt(:), ...
    'VariableNames',{'GateDelay_us','HalfWave','UrOld_V','UrNew_V','ZetaOld','ZetaNew', ...
                     'PlasmaAreaMask_cm2','AreaEq14Old_cm2','AreaEq14New_cm2', ...
                     'MeanExtentOld_mm','MeanExtentNew_mm','TotalInt337'}), ...
           fullfile(OutDir,'Fig13_OldVsNew_Imaging.csv'));
writetable(table(PhaseC,Ub,Vb,Qb,dQdU*1e12,'VariableNames', ...
           {'Phase','Ur_V','V_V','Qp_C','dQpdUr_pF'}),fullfile(OutDir,'Fig13_ChargeDerivative.csv'));
fprintf('Figure and Origin tables written to %s\n',OutDir);

%% Local functions
function [t,V,I] = ReadScope(File,AmpGain,CurrentGain)
% Tektronix MSO64 _ALL.csv: TIME, CH1 (voltage monitor), CH2 (current monitor)
M = readmatrix(File,'NumHeaderLines',11);
t = M(:,1);
V = M(:,2)*AmpGain;
I = M(:,3)*CurrentGain;        % A
end

function f = FitFrequency(t,y,f0)
% Frequency that best fits a sine to y, starting from the nominal f0
Res = @(f) norm(y - [sin(2*pi*f*t) cos(2*pi*f*t) ones(size(t))]*([sin(2*pi*f*t) cos(2*pi*f*t) ones(size(t))]\y));
f = fminsearch(Res,f0,optimset('TolX',1e-3,'Display','off'));
end

function X = Phasor(t,y,f)
% Complex amplitude, cosine reference: y = real(X*exp(1i*2*pi*f*t))
c = [sin(2*pi*f*t) cos(2*pi*f*t) ones(size(t))]\(y - mean(y));
X = c(2) - 1i*c(1);
end

function dy = CircDeriv(y,g)
% Savitzky-Golay first derivative of a periodic signal (as in SDBDCleanVersion.m)
p  = size(g,1);
yy = [y(end-p+1:end); y; y(1:p)];
dd = conv(yy,-g(:,2),'same');
dy = dd(p+1:end-p);
end

function E = MeanLitExtent(SideExt)
% Mean extent of the sides that are lit at each gate (0 if none are)
Lit = SideExt > 0;
E = sum(SideExt.*Lit,1)./max(sum(Lit,1),1);
end
