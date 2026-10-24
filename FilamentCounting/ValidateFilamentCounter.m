%% Validation of the filament counter on injected pulses
%
% Known pulses are added to real captures where nothing is happening (more than 250 us after
% the switch, when the ring-down is far below ignition), and the same detection and
% decomposition as FilimentCounterDuncan.m is run on them. Three kinds of injection:
%   Singles - real isolated pulses from the data (their own shape and charge), one at a time
%   Scan    - template-shaped pulses of 0.1-2 nC, for the detection efficiency against charge
%   Pairs   - two real isolated pulses of the same polarity, 0-320 ns apart, for how close
%             two filaments can be and still be counted separately
% The same records without injection give the false-detection rate in pure noise. Both the
% events (the filament count) and the decomposed pulses (exploratory) are scored.
%
% Inputs:  Output/FilamentCounting/FilamentCounter.mat (run FilimentCounterDuncan.m first)
%          and the 0.12 captures listed in NoiseFiles.
% Outputs: Output/FilamentCounting/Validation (figure, ValidationSummary.csv)
%
% ASSUMPTIONS
%   - A real isolated pulse's "true" charge is its own measured event charge; that carries the
%     noise of one measurement, so the charge scatter shown is an upper bound.
%   - Injected real pulses carry their own recorded noise on top of the record's, so detection
%     here is slightly pessimistic.

clear

%% Section 1 - Options
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
DataDir  = 'E:\MATLAB\MATLAB\Data\2024\VDBDVoltageAndCurrent\SampledData\Sample_1.0mm\AllData';
OutDir   = fullfile(RepoRoot,'Output','FilamentCounting','Validation');
NoiseFiles = "Sample_1.0_0.12_000" + (5:9) + "_Analog.csv";
Spacing  = 8e-6;                                       % s between injected items
NSingles = 150;
QScan    = linspace(0.1,2,200);                         % nC
Seps_ns  = [0 5 10 15 20 30 40 60 80 120 160 240 320];
NPairs   = 15;                                          % per separation
LibQ     = [0.4 3];                                     % nC; library pulses used for injection
Match    = 25e-9;                                       % s; a found pulse this close is a hit
Span     = [-50e-9 800e-9];                             % s; pulses in this window belong to an item
load(fullfile(RepoRoot,'Output','FilamentCounting','FilamentCounter.mat'))
if ~exist(OutDir,'dir'), mkdir(OutDir); end
rng(1)

%% Section 2 - Injection plan
Lib = find(IsoQ >= LibQ(1) & IsoQ <= LibQ(2));
K = [ones(NSingles,1); 2*ones(numel(QScan),1); 3*ones(numel(Seps_ns)*NPairs,1)];
Sep = [zeros(NSingles + numel(QScan),1); repelem(Seps_ns(:),NPairs)*1e-9];
L1 = Lib(randi(numel(Lib),numel(K),1)); L2 = Lib(randi(numel(Lib),numel(K),1));
QTrue = IsoQ(L1);
QTrue(K == 2) = QScan;
QTrue(K == 3) = IsoQ(L1(K == 3)) + IsoQ(L2(K == 3));
Plan = table(K,Sep,L1,L2,QTrue,'VariableNames',{'Kind','Sep','L1','L2','QTrue_nC'});
Plan = Plan(randperm(height(Plan)),:);
for v = ["File" "tInj_us" "Found" "FoundEv" "QFound_nC"], Plan.(v) = nan(height(Plan),1); end
Plan.Hit = false(height(Plan),1);

%% Section 3 - Inject, detect and decompose
Next = 1; FalseEv = 0; FalsePu = 0; QuietTime = 0;
for f = 1:numel(NoiseFiles)
    [t,V,I] = ReadScopeCSV(fullfile(DataDir,NoiseFiles(f)));
    dt = t(2) - t(1);
    V = V*HVRatio; V = V(min(max((1:numel(V)) + round(VoltageLag/dt),1),numel(V)));
    I = I*ProbeGain;
    Ev0 = DetectDischargeEvents(t,V,I,Ccell,struct('KSigma',KSigma));
    i0 = find(t > Ev0.tSwitch + 250e-6,1); i1 = numel(t) - round(30e-6/dt);
    % false detections without injection, in the same stretch
    E0 = Ev0.Events(Ev0.Events.Start >= i0 & Ev0.Events.Start <= i1,:);
    FalseEv = FalseEv + height(E0); QuietTime = QuietTime + (i1 - i0)*dt;
    for e = 1:height(E0)
        FalsePu = FalsePu + height(DecomposeEvent(E0.Pol(e)*(Ev0.Ic(E0.Start(e):E0.Stop(e)) - E0.Base_A(e)), ...
            Tp,Ev0.Sigma,struct('KDet',KDet,'Eta',Eta)));
    end
    Inj = zeros(size(I)); pos = i0 + round(2e-6/dt); Here = [];
    while Next <= height(Plan) && pos + round(Spacing/dt) < i1
        r = Plan(Next,:);
        s = Ev0.Pol(pos); if s == 0, s = 1; end         % real pulses follow dV/dt; so must these
        switch r.Kind
            case 1
                Inj(pos-Pre:pos+Post) = Inj(pos-Pre:pos+Post) + s*IsoSeg(r.L1,:);
            case 2
                Inj = Inj + s*r.QTrue_nC*1e-9*interp1(Tp.tFine,Tp.h,t - t(pos) - rand*dt,'linear',0);
            case 3
                m = round(r.Sep/dt);
                Inj(pos-Pre:pos+Post) = Inj(pos-Pre:pos+Post) + s*IsoSeg(r.L1,:);
                Inj(pos+m-Pre:pos+m+Post) = Inj(pos+m-Pre:pos+m+Post) + s*IsoSeg(r.L2,:);
        end
        Plan.File(Next) = f; Plan.tInj_us(Next) = (t(pos) - Ev0.tSwitch)*1e6;
        Here(end+1) = Next; Next = Next + 1; pos = pos + round(Spacing/dt); %#ok<SAGROW>
    end
    Ev = DetectDischargeEvents(t,V,I + Inj,Ccell,struct('KSigma',KSigma));
    E = Ev.Events;
    Pt = []; Pq = [];
    for e = 1:height(E)
        y = E.Pol(e)*(Ev.Ic(E.Start(e):E.Stop(e)) - E.Base_A(e));
        Comp = DecomposeEvent(y,Tp,Ev.Sigma,struct('KDet',KDet,'Eta',Eta));
        if isempty(Comp), Comp = table((E.Peak(e) - E.Start(e))*dt*1e9,abs(E.Q_nC(e)),'VariableNames',{'tau_ns','q_nC'}); end
        Pt = [Pt; t(E.Start(e)) - Ev.tSwitch + Comp.tau_ns*1e-9]; %#ok<AGROW>
        Pq = [Pq; Comp.q_nC*abs(E.Q_nC(e))/sum(Comp.q_nC)]; %#ok<AGROW>
    end
    Et = t(E.Peak)' - Ev.tSwitch;
    for j = Here
        t0 = Plan.tInj_us(j)*1e-6;
        In = Pt >= t0 + Span(1) & Pt <= t0 + Plan.Sep(j) + Span(2);
        Plan.Found(j) = nnz(In);
        Plan.FoundEv(j) = nnz(Et >= t0 + Span(1) & Et <= t0 + Plan.Sep(j) + Span(2));
        Plan.QFound_nC(j) = sum(Pq(In));
        Plan.Hit(j) = any(abs(Pt - t0) < Match);
    end
    fprintf('%s: %d items injected, %d events found\n',NoiseFiles(f),numel(Here),height(E));
end
Plan = Plan(~isnan(Plan.File),:);

%% Section 4 - Scores
S = Plan(Plan.Kind == 1,:);
fprintf('Real isolated pulses (%d): found %.0f%%, charge found/true median %.2f (IQR %.2f-%.2f)\n', ...
    height(S),100*mean(S.FoundEv >= 1),median(S.QFound_nC./S.QTrue_nC),prctile(S.QFound_nC./S.QTrue_nC,[25 75]));
fprintf('  counted as more than one: events %.0f%%, decomposed pulses %.0f%%\n',100*mean(S.FoundEv > 1),100*mean(S.Found > 1));
fprintf('Pure noise: %d false events and %d false pulses in %.1f ms\n',FalseEv,FalsePu,QuietTime*1e3);
Sc = Plan(Plan.Kind == 2,:);
Edges = 0.1:0.1:2.0; Mid = Edges(1:end-1) + 0.05;
Eff = arrayfun(@(i) mean(Sc.Hit(Sc.QTrue_nC >= Edges(i) & Sc.QTrue_nC < Edges(i+1))),1:numel(Mid));
P = Plan(Plan.Kind == 3,:);
Res = table(Seps_ns(:),'VariableNames',{'Sep_ns'});
for i = 1:numel(Seps_ns)
    r = abs(P.Sep*1e9 - Seps_ns(i)) < 0.5;
    Res.N(i) = nnz(r);
    Res.EventsAs1(i) = mean(P.FoundEv(r) == 1); Res.EventsAs2(i) = mean(P.FoundEv(r) == 2);
    Res.PulsesAs1(i) = mean(P.Found(r) == 1); Res.PulsesAs2(i) = mean(P.Found(r) == 2);
    Res.PulsesAs3Plus(i) = mean(P.Found(r) >= 3);
    Res.QRatio(i) = median(P.QFound_nC(r)./P.QTrue_nC(r));
end
disp(Res)
Q50 = Mid(find(Eff >= 0.5,1)); Q90 = Mid(find(Eff >= 0.9,1));   % first charge bin reaching 50% / 90%
if isempty(Q50), Q50 = NaN; end
if isempty(Q90), Q90 = NaN; end
Summary = table(height(S),mean(S.FoundEv >= 1),mean(S.FoundEv > 1),mean(S.Found > 1),median(S.QFound_nC./S.QTrue_nC), ...
    iqr(S.QFound_nC./S.QTrue_nC),FalseEv/(QuietTime*1e3),FalsePu/(QuietTime*1e3),Q50,Q90, ...
    'VariableNames',{'Singles','FoundFraction','SplitAsEvents','SplitAsPulses','QRatioMedian','QRatioIQR', ...
    'FalseEventsPerMs','FalsePulsesPerMs','Q50pctDetected_nC','Q90pctDetected_nC'});
writetable(Summary,fullfile(OutDir,'ValidationSummary.csv'));
writetable(Res,fullfile(OutDir,'ValidationPairs.csv'));

%% Section 5 - Figure
fig = figure('Color','w','Position',[40 40 1500 450]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile
plot(Mid,Eff,'ko-','MarkerFaceColor','k'); ylim([0 1.05]); box on
xlabel('injected pulse charge / nC'); ylabel('fraction detected')
title('Detection efficiency (template-shaped pulses)')
nexttile
x = max(Res.Sep_ns,2);
plot(x,Res.EventsAs1,'ko-',x,Res.EventsAs2,'ks-','LineWidth',1.3,'MarkerFaceColor','k'); hold on
plot(x,Res.PulsesAs2,'r^--','LineWidth',1); box on
set(gca,'XScale','log'); xlim([1.5 400]); ylim([0 1.05])
xlabel('separation of the two pulses / ns (0 plotted at 2)'); ylabel('fraction of pairs')
legend('events: counted as 1','events: counted as 2','decomposed pulses: counted as 2','Location','west')
title('Two real pulses: how many are counted')
nexttile
histogram(S.QFound_nC./S.QTrue_nC,0:0.05:2,'FaceColor',[.7 .7 .7]); box on
xlabel('found / true charge'); ylabel('isolated real pulses')
title(sprintf('Charge of single pulses: median %.2f',median(S.QFound_nC./S.QTrue_nC)))
SaveFigure(fig,OutDir,'Validation')

function SaveFigure(fig,Folder,Name)
% Save a figure as .png and .fig and leave it open, titled with its file name.
fig.Name = Name;
fig.Visible = 'on';
exportgraphics(fig,fullfile(Folder,[Name '.png']))
savefig(fig,fullfile(Folder,[Name '.fig']))
end
