%% Discharge-location capacitance model of the positive half-wave Qp
%
% Treats every imaged discharge location (each hexagon side, plus the traced
% outer perimeter) as a capacitor that switches on when it ignites. With the
% gap voltage clamped once a location burns, location k adds
%     Qp_k = c * len_k * (Ur - Ur_k)        for Ur > Ur_k
% so the plasma charge derivative is the capacitance of everything lit:
%     dQp/dUr = c * L_lit(Ur),   L_lit = sum of the lit edge lengths
% c (pF per mm of edge) is one number for every location, fitted to the
% measured dQp/dUr of the positive half-wave. The question is whether the
% imaged ignition sequence alone reproduces the structure in that curve.
%
% Models compared (each predicted Qp is passed through exactly the same
% phase-averaged Savitzky-Golay derivative as the measurement):
%   Location         - sides switch on at their imaged ignition gate (1 parameter)
%   Location + burst - as above, plus each location dumps c*len_k*dU at
%                      ignition (overvoltage released at once; 2 parameters).
%                      Steps in L_lit can only make dQp/dUr rise; bursts are
%                      what can make peaks.
%   Lit perimeter    - continuous lit share of every side's edge (zeta-like)
%   Plasma area      - thresholded plasma area (paper Fig 13B), c in pF/mm^2
%
% Timing: the imaging is placed on the scope time axis by the emission-end
% alignment of ICCDFilterRatioEoN.m (light stops at the extremum of the
% electrode voltage Ur). As a check, every model is also refitted over a range
% of extra offsets (OffsetScan). With the earlier alignment on the applied
% voltage, all four models independently preferred +7 us, the lag of Ur behind
% V; with the Ur alignment the best offset should sit near zero.
%
% Inputs: Output/Cache/BigHexDerived.mat (run ICCDFilterRatioEoN.m first) and
% the Big Hex scope traces. Outputs: Output/CapacitanceModel. Figures are left
% open and saved as .png and .fig.

clear

%% Section 1 - Options
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
OutDir   = fullfile(RepoRoot,'Output','CapacitanceModel');
if ~exist(OutDir,'dir'), mkdir(OutDir); end
D = load(fullfile(RepoRoot,'Output','Cache','BigHexDerived.mat'));
ScopeOn  = D.ScopeFile;
ScopeOff = fullfile(fileparts(fileparts(ScopeOn)),'PlasmaOff','bigkex1vpp3000hz_ALL.csv');

PhaseBins        = 1000;   % as in Fig13Compare.m
SGOrder          = 4;
SGWindow         = 41;
TurningPointFrac = 0.10;
CapacitancePerArea = 1.245e-7;   % F/m^2, SDBDCleanVersion.m; only to express c as a width
OffsetScan       = -10:1:16; % us added to the imaging time axis (later = positive)

%% Section 2 - Measured Ur, Qp and dQp/dUr (same pipeline as Fig13Compare.m)
[tOff,VOff,IOff] = ReadScope(ScopeOff,D.AmpGain,D.CurrentGain/1000);
fOff = FitFrequency(tOff,VOff,D.Fdrive);
Z    = Phasor(tOff,VOff,fOff)/Phasor(tOff,IOff,fOff);
Resr = real(Z);
Coff = -1/(2*pi*fOff*imag(Z));

[t,V,I] = ReadScope(ScopeOn,D.AmpGain,D.CurrentGain/1000);
f  = FitFrequency(t,V,D.Fdrive);
P  = 1/f;
Ncyc = floor((t(end) - t(1))/P);
k  = t < t(1) + Ncyc*P;
t = t(k);  V = V(k);  I = I(k) - mean(I(k));
Ur = V - I*Resr;
Q  = cumtrapz(t,I);  Q = Q - mean(Q);
Qp = Q - Coff*Ur;

Phase = mod((t - t(1))/P,1);
b   = min(floor(Phase*PhaseBins) + 1,PhaseBins);
cnt = accumarray(b,1,[PhaseBins 1]);
Ub  = accumarray(b,Ur,[PhaseBins 1])./cnt;
Qb  = accumarray(b,Qp,[PhaseBins 1])./cnt;
[~,g] = sgolay(SGOrder,SGWindow);
dU   = CircDeriv(Ub,g);
Mask = abs(dU) < TurningPointFrac*max(abs(dU));
Meas = CircDeriv(Qb,g)./dU;
Meas(Mask) = NaN;
PhaseC = ((0:PhaseBins-1)' + 0.5)/PhaseBins;

% Rising branch of Ur (positive half-wave), from the Ur minimum to the maximum
[~,iMin] = min(Ub);  [~,iMax] = max(Ub);
Rise = mod((0:PhaseBins-1)' - (iMin-1),PhaseBins) <= mod(iMax - iMin,PhaseBins);
RiseIdx = mod((iMin-1):(iMin-1+mod(iMax-iMin,PhaseBins)),PhaseBins)' + 1;   % in order

% Everything the models need, so they can be refitted for any time offset
Ctx = struct('D',D,'t1',t(1),'P',P,'PhaseBins',PhaseBins,'RiseIdx',RiseIdx,'Rise',Rise, ...
             'Ub',Ub,'dU',dU,'g',g,'Meas',Meas,'Mask',Mask,'PhaseC',PhaseC);
Order = nan(PhaseBins,1);  Order(RiseIdx) = 1:numel(RiseIdx);
Ctx.Order = Order;

%% Section 3 - Fit the models at the imaging alignment and over the offset scan
Mod0 = FitModels(0,Ctx);
Names = {'Location','Burst','Zeta','Area'};
R2 = nan(numel(OffsetScan),numel(Names));
for k = 1:numel(OffsetScan)
    Mk = FitModels(OffsetScan(k),Ctx);
    for m = 1:numel(Names), R2(k,m) = Mk.M.(Names{m}).R2; end
end
[~,kBest] = max(R2,[],1);
BestOffset = OffsetScan(kBest);
ModB = FitModels(BestOffset(1),Ctx);        % at the location model's best offset

for pass = 1:2
    if pass == 1, Mo = Mod0; Lab = sprintf('imaging alignment (emission end on the %s extremum)',D.AlignTo);
    else,         Mo = ModB; Lab = sprintf('best offset for the location model (%+d us)',BestOffset(1)); end
    M = Mo.M;
    fprintf('\nPositive half-wave dQp/dUr, %s:\n',Lab);
    fprintf('  Location model        : c = %.3f pF/mm of edge, R^2 = %.3f (= %.2f mm of dielectric at %.3g F/m^2)\n', ...
            M.Location.Coef(1)*1e12,M.Location.R2,M.Location.Coef(1)*1e3/CapacitancePerArea*1e3,CapacitancePerArea);
    fprintf('  Location + burst      : c = %.3f pF/mm, burst %.0f V per location, R^2 = %.3f\n', ...
            M.Burst.Coef(1)*1e12,M.Burst.Coef(2)/M.Burst.Coef(1),M.Burst.R2);
    fprintf('  Lit perimeter (zeta)  : c = %.3f pF/mm, R^2 = %.3f\n',M.Zeta.Coef(1)*1e12,M.Zeta.R2);
    fprintf('  Plasma area (mask)    : c = %.3f pF/mm^2, R^2 = %.3f\n',M.Area.Coef(1)*1e12,M.Area.R2);
end
fprintf('\nBest extra offset per model (us): location %+d, location+burst %+d, lit perimeter %+d, plasma area %+d\n', ...
        BestOffset);
fprintf('Lit edge length at the end of the half-wave: %.0f mm (%d of %d sides ignited)\n', ...
        max(Mod0.LLoc),nnz(Mod0.Ignited),numel(Mod0.Ignited));

%% Section 4 - Figure
fig = figure('Color','w','Position',[40 40 1500 950]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
x = Ub(RiseIdx);
Mods = {Mod0,ModB};
Ttl  = {'Imaging alignment (as used everywhere else)', ...
        sprintf('Imaging shifted %+d \\mus (best fit of the location model)',BestOffset(1))};
for c = 1:2
    M = Mods{c}.M;
    nexttile(c,[2 1]); hold on
    plot(x,Meas(RiseIdx)*1e12,'k-','LineWidth',2.2)
    plot(x,M.Location.Pred(RiseIdx)*1e12,'-','Color',[0.15 0.35 0.85],'LineWidth',1.6)
    plot(x,M.Burst.Pred(RiseIdx)*1e12,'-','Color',[0.85 0.25 0.1],'LineWidth',1.3)
    plot(x,M.Zeta.Pred(RiseIdx)*1e12,'--','Color',[0.2 0.6 0.3],'LineWidth',1.2)
    plot(x,M.Area.Pred(RiseIdx)*1e12,':','Color',[0.4 0.4 0.4],'LineWidth',1.6)
    ylabel('dQ_p/dU_r / pF'); box on; grid on; ylim([-50 600])
    legend({'Measured', ...
            sprintf('Location, fixed C per side (c = %.2f pF/mm, R^2 = %.2f)',M.Location.Coef(1)*1e12,M.Location.R2), ...
            sprintf('Location + ignition burst (%.0f V, R^2 = %.2f)',M.Burst.Coef(2)/M.Burst.Coef(1),M.Burst.R2), ...
            sprintf('Lit perimeter (R^2 = %.2f)',M.Zeta.R2), ...
            sprintf('Plasma area (R^2 = %.2f)',M.Area.R2)},'Location','northwest','FontSize',8)
    title(Ttl{c})
    xl = xlim;
    nexttile(4 + c); hold on
    bar(Mods{c}.GateUr,Mods{c}.NIgn,1,'FaceColor',[0.85 0.25 0.1],'EdgeColor','none','FaceAlpha',0.6)
    ylabel('Sides igniting per gate')
    yyaxis right
    plot(x,Mods{c}.LLoc(RiseIdx),'-','Color',[0.15 0.35 0.85],'LineWidth',1.5)
    ylabel('Lit edge length / mm')
    xlabel('U_r / V'); xlim(xl); box on
end
SaveFigure(fig,OutDir,'LocationCapacitanceModel')

fig2 = figure('Color','w','Position',[80 80 700 450]);
plot(OffsetScan,R2,'o-','LineWidth',1.4); grid on; box on
xlabel('Extra offset added to the imaging time axis / \mus'); ylabel('R^2 of the fit to dQ_p/dU_r')
legend({'Location','Location + burst','Lit perimeter','Plasma area'},'Location','best')
title('How well each model fits, against the imaging-to-scope timing')
SaveFigure(fig2,OutDir,'LocationCapacitanceModel_OffsetScan')

writetable(array2table([OffsetScan(:) R2],'VariableNames',[{'ExtraOffset_us'} strcat('R2_',Names)]), ...
           fullfile(OutDir,'LocationCapacitanceModel_OffsetScan.csv'));
Out = table;
for c = 1:2
    M = Mods{c}.M;
    Offs = [0 BestOffset(1)];
    Out = [Out; table(repmat(Offs(c),4,1), ...
        {'Location';'Location + burst';'Lit perimeter';'Plasma area'}, ...
        [M.Location.Coef(1);M.Burst.Coef(1);M.Zeta.Coef(1);M.Area.Coef(1)]*1e12, ...
        [NaN;M.Burst.Coef(2)/M.Burst.Coef(1);NaN;NaN], ...
        [M.Location.R2;M.Burst.R2;M.Zeta.R2;M.Area.R2], ...
        'VariableNames',{'ExtraOffset_us','Model','c_pF_per_mm_or_mm2','Burst_V','R2'})]; %#ok<AGROW>
end
writetable(Out,fullfile(OutDir,'LocationCapacitanceModel_Fits.csv'));
writetable(table(PhaseC(RiseIdx),Ub(RiseIdx),Meas(RiseIdx)*1e12, ...
    Mod0.M.Location.Pred(RiseIdx)*1e12,ModB.M.Location.Pred(RiseIdx)*1e12, ...
    Mod0.M.Area.Pred(RiseIdx)*1e12,ModB.M.Area.Pred(RiseIdx)*1e12,Mod0.LLoc(RiseIdx),ModB.LLoc(RiseIdx), ...
    'VariableNames',{'Phase','Ur_V','Measured_pF','Location_pF','Location_BestOffset_pF', ...
    'PlasmaArea_pF','PlasmaArea_BestOffset_pF','LitLength_mm','LitLength_BestOffset_mm'}), ...
    fullfile(OutDir,'LocationCapacitanceModel.csv'));

%% Local functions
function Mod = FitModels(ExtraUs,C)
% Builds the lit length of every location on the rising branch for a given
% extra offset of the imaging time axis, and fits all four models.
D = C.D;  RiseIdx = C.RiseIdx;  Order = C.Order;  PB = C.PhaseBins;
GatePhase = mod(((D.GateDelay + D.GateWidth/2)*1e-6 + D.ScopeTimeShift + ExtraUs*1e-6 - C.t1)/C.P,1);
gPos = [];
for h = 1:numel(D.HumpStart)
    gg = D.HumpStart(h):D.HumpEnd(h);
    if D.HalfWave(gg(1)) > 0, gPos = gg; end
end
OrderAt = @(ph) Order(min(floor(ph*PB) + 1,PB));

% Sides switch fully on at their imaged ignition gate and stay on until the
% Ur maximum (the extent does not recede, Fig 10)
Tp = D.SideTable(strcmp(D.SideTable.HalfWave,'Positive'),:);
IgnGate = round(Tp.IgnitionTime_us/D.GateWidth) + 1;
Ignited = ~isnan(IgnGate);
IgnOrder = nan(size(IgnGate));
IgnOrder(Ignited) = OrderAt(GatePhase(IgnGate(Ignited)));
len = D.SideLenMM(:);

GateOrd   = OrderAt(GatePhase(gPos));
CarryRise = @(GateVals) CarryOnRise(GateVals,GateOrd,Order,RiseIdx);
ArcOnBins = CarryRise(D.ArcLit(gPos));        % outer perimeter, continuous lit share

LLoc = zeros(PB,1);  LBurst = zeros(PB,1);
for j = RiseIdx'
    on = Ignited & IgnOrder <= Order(j);
    LBurst(j) = sum(len(on));
    LLoc(j)   = LBurst(j) + D.ArcLenMM*ArcOnBins(j);
end
LitLen = sum(D.SideLit(:,gPos).*len,1,'omitnan') + D.ArcLenMM*D.ArcLit(gPos);   % sides off the image have no edge
LZeta = CarryRise(LitLen);
LArea = CarryRise(D.PlasmaArea(gPos));

% The model charge (per unit c) is built on the rising branch, closed back to
% zero along the falling branch so the cycle has no jump, and differentiated
% with the same circular Savitzky-Golay filter and turning-point mask as the
% measurement. Every model is linear in its coefficients.
Use = C.Rise & ~isnan(C.Meas);
Fit = @(Basis) LinearFit(C.Meas,Basis,Use,C.Mask);
BLoc = DerivOf(Integrate(LLoc,C.Ub,RiseIdx),C.Ub,RiseIdx,C.dU,C.g);
Mod.M.Location = Fit(BLoc);
Mod.M.Burst    = Fit([BLoc DerivOf(Step(LBurst,RiseIdx),C.Ub,RiseIdx,C.dU,C.g)]);
Mod.M.Zeta     = Fit(DerivOf(Integrate(LZeta,C.Ub,RiseIdx),C.Ub,RiseIdx,C.dU,C.g));
Mod.M.Area     = Fit(DerivOf(Integrate(LArea,C.Ub,RiseIdx),C.Ub,RiseIdx,C.dU,C.g));
Mod.LLoc    = LLoc;
Mod.Ignited = Ignited;
Mod.GateUr  = interp1([C.PhaseC-1; C.PhaseC; C.PhaseC+1],[C.Ub; C.Ub; C.Ub],GatePhase(gPos));
Mod.NIgn    = arrayfun(@(gi) nnz(IgnGate == gi),gPos);
end

function Lb = CarryOnRise(Vals,GateOrd,Order,RiseIdx)
% Value of the latest gate at or before each rising-branch bin
Lb = zeros(numel(Order),1);
ok = ~isnan(GateOrd(:));
Vals = Vals(:);
[go,ix] = sort(GateOrd(ok));
v = Vals(ok);  v = v(ix);
for j = RiseIdx'
    k = find(go <= Order(j),1,'last');
    if ~isempty(k), Lb(j) = v(k); end
end
end

function Qm = Integrate(L,Ub,RiseIdx)
% Charge added along the rising branch for unit c: cumulative sum of L*dUr
Qm = zeros(size(Ub));
dUr = diff(Ub(RiseIdx));
Qm(RiseIdx) = [0; cumsum(L(RiseIdx(2:end)).*dUr)];
end

function Qm = Step(L,RiseIdx)
% Charge released at once at ignition, for unit c*dU
Qm = zeros(size(L));
Qm(RiseIdx) = L(RiseIdx);
end

function y = DerivOf(Qrise,Ub,RiseIdx,dU,g)
% dQ/dUr of a model charge defined on the rising branch (0 at the Ur minimum).
% On the falling branch the charge is taken back to 0 in proportion to Ur, so
% the periodic derivative has no jump; that branch is not compared.
Qm = zeros(size(Ub));
Qm(RiseIdx) = Qrise(RiseIdx);
Fall = setdiff((1:numel(Ub))',RiseIdx);
Umin = Ub(RiseIdx(1));  Umax = Ub(RiseIdx(end));
Qm(Fall) = Qrise(RiseIdx(end))*(Ub(Fall) - Umin)/(Umax - Umin);
y = CircDeriv(Qm,g)./dU;
end

function F = LinearFit(y,X,Use,Mask)
% y ~ X*Coef, least squares over Use; predictions blanked where the
% measurement is (turning points)
F.Coef = X(Use,:)\y(Use);
F.Pred = X*F.Coef;
F.Pred(Mask) = NaN;
Res = y(Use) - F.Pred(Use);
F.R2 = 1 - sum(Res.^2)/sum((y(Use) - mean(y(Use))).^2);
end

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

function SaveFigure(fig,Folder,Name)
% Save a figure as .png and .fig and leave it open, titled with its file name.
fig.Name = Name;
fig.Visible = 'on';
exportgraphics(fig,fullfile(Folder,[Name '.png']))
savefig(fig,fullfile(Folder,[Name '.fig']))
end
