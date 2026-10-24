%% Filament count and charge per filament from the VDBD ring-down captures
%
% Each capture (Tektronix MSO64, 1.6 ns sampling, 2 ms record) holds one excitation of the
% supply: a ~27 kHz damped ring-down. Discharges occur only in the first ~110 us, in every
% half-cycle while the amplitude is high enough, with the polarity of dV/dt. CH1 is the HV
% probe (HVRatio), CH2 the Pearson current monitor (ProbeGain).
%
% Steps
%   2. Cell capacitance check. Duty 0.08 never ignites, so its current is pure displacement
%      current: Ccell and the time lag between the two probes are fitted there, by voltage
%      window (a Ccell that grew with amplitude would mean hidden discharge current). The
%      same fit on the discharge-free stretches of every other capture is the check.
%   3. Events. Ccell*dV/dt and the CH2 offset are removed and pulses above KSigma times the
%      noise are found and integrated against a local baseline (DetectDischargeEvents).
%   4. Pulses (exploratory). Each event is fitted with shifted copies of the measured
%      single-pulse response (BuildPulseTemplate, DecomposeEvent). On injected real pulses
%      (ValidateFilamentCounter.m) this splits a quarter of single pulses and miscounts pairs
%      at every separation tried (0-320 ns) and every sample interval (1.6-12.8 ns): the pulse
%      shape varies too much (17%) and big events ring differently. Pulse counts are reported
%      but not used as the filament count.
%   5. Filaments. Event charges are tested for humps at multiples of a single-filament charge
%      q1 (FitChargeQuanta). Only if the test supports it does an event of charge Q count as
%      round(Q/q1) filaments; otherwise each event counts as one (a lower bound: filaments less
%      than ~100-300 ns apart, or too small to detect, are not separated).
%   6. Half-cycles and Lissajous. Each ring-down is split into half-cycles at its voltage
%      peaks. The integrated current Q (CH2 minus its mean over the whole pulse) drifts with
%      the slow CH2 baseline wander, so it is corrected on the stretch after each voltage peak,
%      where no discharge can start and Q - Ccell*V must stay flat; the step in that level
%      across a half-cycle is the half-cycle's total discharge charge, including filaments too
%      small to detect. Comparing it with the sum of event charges gives that hidden share.
%   7. Dielectric capacitance and gap voltage (equivalent circuit, Pipa & Brandenburg).
%      zeta_diel, the Q-V slope while the discharge is on, comes from the Lissajous charge:
%      zeta = Ccell + dQ_discharge/dV from the first event to the voltage peak. With only part
%      of the electrode discharged, zeta_diel <= C_d. The true C_d needs the dielectric
%      (EpsR, DielThick_mm); if given, the gap voltage at every event U_g = kappa*V - Q_d/C_d,
%      the gap charge Q/kappa and the energy per event follow, kappa = 1 - Ccell/C_d.
%   8. The ring-down decays, so every capture is a voltage sweep: results are tabulated per
%      half-cycle against its peak voltage, and per capture and setting.
%
% Inputs:  E:\MATLAB\MATLAB\Data\2024\VDBDVoltageAndCurrent\SampledData\Sample_1.0mm\AllData
%          (outside the repo)
% Outputs: Output/FilamentCounting: figures (.png, .fig), tables (.csv) and FilamentCounter.mat
%          (calibration, template, quanta fit) for ValidateFilamentCounter.m
%
% ASSUMPTIONS
%   - CH2 is 1 A/V (ProbeGain) and CH1 1000:1 (HVRatio). Pearson monitors are specified into
%     1 MOhm; into 50 Ohm their sensitivity halves, which would scale every charge.
%   - The Average-mode folders are not used: averaging erases the random-in-time pulses.
%   - Filaments too small to clear the threshold, and filaments merged into one event, are
%     not counted separately: counts are lower bounds. The small ones' charge is in the
%     Lissajous charge but not in the event charges.
%   - Gap voltages assume no surface charge is left from the previous excitation, and are
%     electrode averages: right after a filament the gap voltage at its own spot is lower.

clear

%% Section 1 - Options
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
DataDir  = 'E:\MATLAB\MATLAB\Data\2024\VDBDVoltageAndCurrent\SampledData\Sample_1.0mm\AllData';
OutDir   = fullfile(RepoRoot,'Output','FilamentCounting');
Pattern     = 'Sample_1.0_%s_*_Analog.csv';
Settings    = {'0.08','0.09','0.10','0.11','0.12'};   % supply duty settings
CalSettings = {'0.08'};                                % settings that never ignite
HVRatio   = 1000;              % CH1: V per probe V
ProbeGain = 1;                 % CH2: A per V
EpsR      = NaN;               % dielectric relative permittivity   } leave NaN if unknown:
DielThick_mm = NaN;            % dielectric thickness (mm)            } gap voltage is skipped
Gap_mm    = 1.0;               % gas gap (mm)
CalWin_us = [3 20; 20 60; 60 150; 150 400];            % capacitance fit windows after the switch
CalUse    = [2 3];             % windows that set Ccell: strong signal, past the first swing
KSigma    = 6;                 % event threshold, x noise sigma
KDet      = 6;                 % pulse detection in the decomposition, matched-filter sigma
Eta       = 0.5;               % a new pulse must exceed this share of the fit around it
TemplateQ = [0.3 1.0];         % nC; isolated events in this range build the pulse template
IsoTime   = 1e-6;              % s; an isolated event has no other event this close
Pre = 25; Post = 500;          % template window, samples (-40 ns .. +800 ns)
QRange    = [0.7 6];           % nC; charge-quanta fit range (detection is incomplete below 0.7 nC)
NBoot     = 200;               % bootstrap replicates of the charge-quanta test
PQuanta   = 0.05;              % quantization is used for counting only if p_boot < PQuanta
NBootZ    = 1000;              % bootstrap replicates (over captures) for zeta_diel
WaveStep  = 4;                 % keep every 4th sample of the ring-down for Sections 6-7
WaveSpan  = [-5e-6 400e-6];    % s after the switch
RepRate   = NaN;               % excitation repetition rate (Hz), if known: adds rates per second
ExampleFile = 'Sample_1.0_0.12_0000_Analog.csv';
if ~exist(OutDir,'dir'), mkdir(OutDir); end

%% Section 2 - Cell capacitance check
% No discharge: I(t) = Ccell*dV/dt(t + VoltageLag), i.e. the voltage channel reads VoltageLag late.
Cal = ListCaptures(DataDir,Pattern,CalSettings);
LagScan = (-250:5:250)*1e-9;
NW = size(CalWin_us,1);
CalRms = zeros(height(Cal),numel(LagScan));
CalCS  = zeros(height(Cal),NW,numel(LagScan));
CalV   = zeros(height(Cal),NW);
for f = 1:height(Cal)
    [t,V,I] = ReadScopeCSV(Cal.File{f});
    V = V*HVRatio; I = I*ProbeGain; dt = t(2)-t(1);
    dVdt = SmoothDeriv(V,dt);
    iSw = SwitchIndex(t,V);
    I = I - mean(I);                                   % DC bias: mean over the whole pulse
    m = find(t > t(iSw) + CalWin_us(CalUse(1),1)*1e-6 & t < t(iSw) + CalWin_us(CalUse(end),2)*1e-6);
    mw = cell(1,NW);
    for w = 1:NW
        mw{w} = find(t > t(iSw) + CalWin_us(w,1)*1e-6 & t < t(iSw) + CalWin_us(w,2)*1e-6);
        CalV(f,w) = max(abs(V(mw{w})));
    end
    for j = 1:numel(LagScan)
        k = round(LagScan(j)/dt);
        c = dVdt(m+k)' \ I(m)';
        CalRms(f,j) = rms(I(m) - c*dVdt(m+k));
        for w = 1:NW, CalCS(f,w,j) = dVdt(mw{w}+k)' \ I(mw{w})'; end
    end
end
[~,jb] = min(mean(CalRms./min(CalRms,[],2),1));
VoltageLag = LagScan(jb);
CalC = CalCS(:,:,jb);
Ccell = mean(CalC(:,CalUse),'all');
CalTable = table(Cal.Name,repmat(VoltageLag*1e9,height(Cal),1),CalC*1e12,CalV/1e3,squeeze(CalCS(:,:,LagScan == 0))*1e12, ...
    'VariableNames',{'File','VoltageLag_ns','C_pF','Vmax_kV','CNoLag_pF'});
writetable(splitvars(CalTable),fullfile(OutDir,'CellCapacitanceCheck.csv'));
fprintf('Cell capacitance check (%d discharge-free captures): voltage channel lags by %.0f ns\n',height(Cal),VoltageLag*1e9);
for w = 1:NW
    fprintf('  %3.0f-%3.0f us after switch, |V| up to %.1f kV: C = %.3f +- %.3f pF\n',CalWin_us(w,:), ...
        mean(CalV(:,w))/1e3,mean(CalC(:,w))*1e12,std(CalC(:,w))*1e12);
end
fprintf('  Ccell = %.3f pF (the old script used 7.799 pF)\n',Ccell*1e12);

%% Section 3 - Discharge events in every capture
Caps = ListCaptures(DataDir,Pattern,Settings);
NC = height(Caps);
for v = ["VScale" "Sigma_mA" "Offset_mA" "PulseMean_mA" "CQuiet_pF"], Caps.(v) = nan(NC,1); end
Events = table(); Segs = {}; IsoSeg = zeros(0,Pre+Post+1); IsoQ = zeros(0,1);
Wave = cell(NC,1);
for f = 1:NC
    [t,V,I,Info] = ReadScopeCSV(Caps.File{f});
    dt = t(2)-t(1);
    V = ShiftSamples(V*HVRatio,-round(VoltageLag/dt));
    I = I*ProbeGain;
    Ev = DetectDischargeEvents(t,V,I,Ccell,struct('KSigma',KSigma));
    E = Ev.Events;
    Caps.VScale(f) = Info.VScale(2); Caps.Sigma_mA(f) = Ev.Sigma*1e3;
    Caps.Offset_mA(f) = Ev.Offset*1e3; Caps.PulseMean_mA(f) = mean(I)*1e3;
    E.Capture = repmat(f,height(E),1);
    E.tStart_us = (t(E.Start)' - Ev.tSwitch)*1e6;
    for e = 1:height(E)
        Segs{end+1,1} = E.Pol(e)*(Ev.Ic(E.Start(e):E.Stop(e)) - E.Base_A(e)); %#ok<SAGROW>
        Others = E.Peak([1:e-1 e+1:end]);
        if all(abs(Others - E.Peak(e))*dt > IsoTime) && E.Peak(e) > Pre + 200 && E.Peak(e) + Post <= numel(t)
            bw = E.Start(e) - round(200e-9/dt):E.Start(e) - 1;
            IsoSeg(end+1,:) = E.Pol(e)*(Ev.Ic(E.Peak(e)-Pre:E.Peak(e)+Post) - mean(Ev.Ic(bw))); %#ok<SAGROW>
            IsoQ(end+1,1) = abs(E.Q_nC(e)); %#ok<SAGROW>
        end
    end
    Events = [Events; E]; %#ok<AGROW>
    % Capacitance check on this capture's own discharge-free stretches (3-150 us, +-2 us
    % around events removed). Sub-threshold discharge current would push it above Ccell.
    Busy = false(size(t));
    for e = 1:height(E), Busy(max(1,E.Start(e)-round(2e-6/dt)):min(numel(t),E.Stop(e)+round(2e-6/dt))) = true; end
    m = find(t > Ev.tSwitch + 3e-6 & t < Ev.tSwitch + 150e-6 & ~Busy & ~Ev.Masked);
    Caps.CQuiet_pF(f) = Ev.dVdt(m)' \ (I(m) - mean(I))'*1e12;
    % Ring-down kept for Sections 6-7: V and the integrated current (DC bias = mean over the pulse)
    w = find(t - Ev.tSwitch >= WaveSpan(1) & t - Ev.tSwitch <= WaveSpan(2));
    Q = cumtrapz(t(w),I(w) - mean(I));
    k = 1:WaveStep:numel(w);
    Wave{f} = struct('t',t(w(k)) - Ev.tSwitch,'V',V(w(k)),'Q',Q(k) - Q(1));
    fprintf('%s: %d events, noise %.2f mA, CH2 offset %.2f mA (pulse mean %.2f mA)\n',Caps.Name{f}, ...
        height(E),Caps.Sigma_mA(f),Caps.Offset_mA(f),Caps.PulseMean_mA(f));
end

%% Section 4 - Pulse template and decomposition of every event
Use = IsoQ >= TemplateQ(1) & IsoQ <= TemplateQ(2);
Tp = BuildPulseTemplate(IsoSeg(Use,:),dt);
Half = find(Tp.h > max(Tp.h)/2);
TpFWHM = (Half(end) - Half(1))*dt/Tp.Up;
fprintf('Pulse template: %d isolated pulses of %.1f-%.1f nC, FWHM %.1f ns, shape spread %.0f%%\n', ...
    Tp.N,TemplateQ,TpFWHM*1e9,100*Tp.Spread);
PCols = {'Capture','Event','Pol','t_us','q_nC','qFit_nC','V_kV'};
Pulses = cell(height(Events),1);
for e = 1:height(Events)
    f = Events.Capture(e);
    Comp = DecomposeEvent(Segs{e},Tp,Caps.Sigma_mA(f)/1e3,struct('KDet',KDet,'Eta',Eta));
    if isempty(Comp)                                   % too small to fit: keep it as one pulse
        Comp = table((Events.Peak(e) - Events.Start(e))*dt*1e9,abs(Events.Q_nC(e)),'VariableNames',{'tau_ns','q_nC'});
    end
    n = height(Comp);
    tp = Events.tStart_us(e) + Comp.tau_ns/1e3;
    qs = Events.Pol(e)*Comp.q_nC*abs(Events.Q_nC(e))/sum(Comp.q_nC);
    Pulses{e} = table(repmat(f,n,1),repmat(e,n,1),repmat(Events.Pol(e),n,1),tp,qs,Events.Pol(e)*Comp.q_nC, ...
        interp1(Wave{f}.t,Wave{f}.V,tp*1e-6)/1e3,'VariableNames',PCols);
end
Pulses = vertcat(Pulses{:});
if isempty(Pulses), Pulses = cell2table(cell(0,numel(PCols)),'VariableNames',PCols); end
Pulses.Setting = Caps.Setting(Pulses.Capture);
Events.Setting = Caps.Setting(Events.Capture);
EvFit = accumarray(Pulses.Event,abs(Pulses.qFit_nC),[height(Events) 1]);
fprintf('Decomposition: %d events -> %d pulses; fitted / measured event charge: median %.2f (IQR %.2f-%.2f)\n', ...
    height(Events),height(Pulses),median(EvFit./abs(Events.Q_nC)),prctile(EvFit./abs(Events.Q_nC),[25 75]));

%% Section 5 - Filament multiplicity from charge quanta
Fq = FitChargeQuanta(abs(Events.Q_nC),QRange,8,NBoot);
UseQuanta = Fq.PBoot < PQuanta && Fq.BICQuanta < Fq.BICSmooth;
if UseQuanta
    Events.N = max(1,round(abs(Events.Q_nC)/Fq.q1)); CountRule = sprintf('multiplicity from q_1 = %.2f nC',Fq.q1);
else
    Events.N = ones(height(Events),1); CountRule = 'no charge quanta: each event counted once';
end
fprintf('Charge quanta (%d events in %.1f-%.1f nC): best q1 = %.2f nC; BIC smooth - quantized = %.1f; bootstrap p = %.3f -> %s\n', ...
    Fq.NFit,QRange,Fq.q1,Fq.BICSmooth - Fq.BICQuanta,Fq.PBoot,CountRule);

%% Section 6 - Half-cycles and Lissajous
HC = cell(NC,1); Rc = cell(NC,1);
for f = 1:NC
    E = Events(Events.Capture == f,:);
    tEv = [E.tStart_us, E.tStart_us + (E.Stop - E.Start)*dt*1e6]*1e-6;
    [H,Rc{f}] = SplitHalfCycles(Wave{f},Ccell,tEv);
    P = Pulses(Pulses.Capture == f,:);
    for v = ["Events" "Pulses" "Filaments"], H.(v) = zeros(height(H),1); end
    for v = ["QApp_nC" "tOnset_us" "VOnset_kV" "QBefore_nC"], H.(v) = nan(height(H),1); end
    for k = 1:height(H)
        inE = E.tStart_us >= H.tWin_us(k,1) & E.tStart_us < H.tWin_us(k,2);
        inP = P.t_us >= H.tWin_us(k,1) & P.t_us < H.tWin_us(k,2);
        H.Events(k) = nnz(inE); H.Pulses(k) = nnz(inP); H.Filaments(k) = sum(E.N(inE));
        H.QApp_nC(k) = sum(E.Q_nC(inE));
        % onset: first event of the half-cycle's own polarity, and the charge moved before it
        j = find(inE & E.Pol == H.Pol(k),1);
        if ~isempty(j)
            H.tOnset_us(k) = E.tStart_us(j); H.VOnset_kV(k) = E.V_kV(j);
            H.QBefore_nC(k) = sum(E.Q_nC(E.Start < E.Start(j)));
        end
    end
    H.Capture = repmat(f,height(H),1); H.Setting = repmat(Caps.Setting(f),height(H),1);
    H.Shot = repmat(Caps.Shot(f),height(H),1);
    HC{f} = H;
end
HC = vertcat(HC{:});
HC.QPerFil_nC = abs(HC.QApp_nC)./HC.Filaments;
Quiet = HC.Events == 0 & HC.Setting > min(cellfun(@str2double,CalSettings)) & ~isnan(HC.QMeas_nC);
Busy = HC.Events > 0 & ~isnan(HC.QMeas_nC);
Hidden = HC.Pol(Busy).*(HC.QMeas_nC(Busy) - HC.QApp_nC(Busy));
fprintf('Lissajous check: half-cycles without events %.2f +- %.2f nC (n = %d, scatter %.1f nC per shot)\n', ...
    mean(HC.Pol(Quiet).*HC.QMeas_nC(Quiet)),std(HC.QMeas_nC(Quiet))/sqrt(nnz(Quiet)),nnz(Quiet),std(HC.QMeas_nC(Quiet)));
fprintf('  half-cycles with events: Lissajous - event charge = %.2f +- %.2f nC; Lissajous / event charge = %.2f\n', ...
    mean(Hidden),std(Hidden)/sqrt(numel(Hidden)),sum(HC.Pol(Busy).*HC.QMeas_nC(Busy))/sum(HC.Pol(Busy).*HC.QApp_nC(Busy)));

%% Section 7 - Dielectric capacitance and gap voltage
% zeta_diel = Ccell + (discharge charge)/(voltage swing from the first pulse to the peak),
% as a ratio of sums over the half-cycles of a setting (one half-cycle is far too noisy),
% with a bootstrap over captures. Both the Lissajous charge and the event charge are used.
SetVal = cellfun(@str2double,Settings);
Zt = table();
for s = [SetVal NaN]
    if isnan(s), r = Busy & ~isnan(HC.VOnset_kV); else, r = Busy & ~isnan(HC.VOnset_kV) & HC.Setting == s; end
    if nnz(r) < 3, continue, end
    Hr = HC(r,:);
    Swing = Hr.Pol.*(Hr.Vpk_kV - Hr.VOnset_kV)*1e3;
    Zl = Ccell + sum(Hr.Pol.*Hr.QMeas_nC)*1e-9/sum(Swing);
    Zp = Ccell + sum(Hr.Pol.*Hr.QApp_nC)*1e-9/sum(Swing);
    Cp = unique(Hr.Capture); ZB = nan(NBootZ,1);
    for b = 1:NBootZ
        pick = Cp(randi(numel(Cp),numel(Cp),1));
        i = cell2mat(arrayfun(@(c) find(Hr.Capture == c),pick,'UniformOutput',false));
        ZB(b) = Ccell + sum(Hr.Pol(i).*Hr.QMeas_nC(i))*1e-9/sum(Swing(i));
    end
    Zt = [Zt; table(s,nnz(r),Zl*1e12,prctile(ZB,2.5)*1e12,prctile(ZB,97.5)*1e12,Zp*1e12, ...
        'VariableNames',{'Setting','HalfCycles','ZetaLissajous_pF','ZetaLo_pF','ZetaHi_pF','ZetaEvents_pF'})]; %#ok<AGROW>
end
disp(Zt)
writetable(Zt,fullfile(OutDir,'DielectricCapacitance.csv'));
% True C_d from the dielectric, if known: kappa = C_d/(C_d + C_g) per unit area
Kappa = (EpsR/DielThick_mm)/(EpsR/DielThick_mm + 1/Gap_mm);
Cd = Ccell/(1 - Kappa);
GapOK = ~isnan(Kappa);
for v = ["Ug_kV" "QGap_nC" "E_uJ"], Events.(v) = nan(height(Events),1); end
if GapOK
    fprintf('Dielectric: kappa = %.3f, C_d = %.1f pF, C_g = %.1f pF\n',Kappa,Cd*1e12,Ccell*Cd/(Cd - Ccell)*1e12);
    for f = 1:NC
        j = find(Events.Capture == f);                 % events are in time order
        Before = [0; cumsum(Events.Q_nC(j(1:end-1)))]*1e-9;
        Ug = Kappa*Events.V_kV(j)*1e3 - Before/Cd;
        Events.Ug_kV(j) = Ug/1e3;
        Events.QGap_nC(j) = Events.Q_nC(j)/Kappa;
        Events.E_uJ(j) = abs(Events.QGap_nC(j))*1e-9.*max(abs(Ug) - abs(Events.Q_nC(j))*1e-9/Cd/2,0)*1e6;
    end
else
    fprintf('Dielectric permittivity/thickness not set: gap voltage, gap charge and energy are skipped\n');
end
HC.UgOnset_kV = (Kappa*HC.VOnset_kV*1e3 - HC.QBefore_nC*1e-9/Cd)/1e3;
HC.QGap_nC = zeros(height(HC),1); HC.Energy_uJ = zeros(height(HC),1);
for k = 1:height(HC)
    inE = Events.Capture == HC.Capture(k) & Events.tStart_us >= HC.tWin_us(k,1) & Events.tStart_us < HC.tWin_us(k,2);
    HC.QGap_nC(k) = sum(Events.QGap_nC(inE)); HC.Energy_uJ(k) = sum(Events.E_uJ(inE));
end

%% Section 8 - Tables
for v = ["Events" "PulsesExploratory" "Filaments" "PosFilaments" "NegFilaments" "QEvents_nC" "QLissajous_nC" "QGap_nC" "Energy_uJ"]
    Caps.(v) = zeros(NC,1);
end
for f = 1:NC
    E = Events(Events.Capture == f,:);
    Caps.Events(f) = height(E);
    Caps.PulsesExploratory(f) = nnz(Pulses.Capture == f);
    Caps.Filaments(f) = sum(E.N);
    Caps.PosFilaments(f) = sum(E.N(E.Pol > 0)); Caps.NegFilaments(f) = sum(E.N(E.Pol < 0));
    Caps.QEvents_nC(f) = sum(abs(E.Q_nC));
    r = HC.Capture == f & HC.Events > 0;
    Caps.QLissajous_nC(f) = sum(HC.Pol(r).*HC.QMeas_nC(r),'omitnan');
    Caps.QGap_nC(f) = sum(abs(E.QGap_nC)); Caps.Energy_uJ(f) = sum(E.E_uJ);
end
writetable(removevars(Caps,'File'),fullfile(OutDir,'FilamentCountPerCapture.csv'));
writetable(removevars(Events,{'Start','Stop','Peak'}),fullfile(OutDir,'Events.csv'));
writetable(Pulses,fullfile(OutDir,'PulsesExploratory.csv'));
writetable(splitvars(HC),fullfile(OutDir,'HalfCycles.csv'));
St = table();
for s = 1:numel(Settings)
    c = Caps.Setting == SetVal(s);
    E = Events(ismember(Events.Capture,find(c)),:);
    Row = table(string(Settings{s}),nnz(c),mean(Caps.Filaments(c)),std(Caps.Filaments(c)), ...
        mean(Caps.PosFilaments(c)),mean(Caps.NegFilaments(c)),mean(Caps.PulsesExploratory(c)), ...
        sum(abs(E.Q_nC))/max(sum(E.N),1),median(abs(E.Q_nC)./E.N),prctile(abs(E.Q_nC)./E.N,25),prctile(abs(E.Q_nC)./E.N,75), ...
        mean(Caps.QEvents_nC(c)),mean(Caps.QLissajous_nC(c)),std(Caps.QLissajous_nC(c))/sqrt(nnz(c)), ...
        sum(abs(E.QGap_nC))/max(sum(E.N),1),mean(Caps.Energy_uJ(c)),mean(Caps.CQuiet_pF(c)),std(Caps.CQuiet_pF(c)), ...
        'VariableNames',{'Setting','Captures','Filaments','FilamentsSD','PosFilaments','NegFilaments','PulsesExploratory', ...
        'QPerFilamentMean_nC','QPerFilamentMedian_nC','QPerFilamentP25_nC','QPerFilamentP75_nC', ...
        'QEventsPerCapture_nC','QLissajousPerCapture_nC','QLissajousSE_nC', ...
        'QGapPerFilament_nC','EnergyPerCapture_uJ','CQuiet_pF','CQuietSD_pF'});
    if ~isnan(RepRate), Row.FilamentsPerSecond = Row.Filaments*RepRate; end
    St = [St; Row]; %#ok<AGROW>
end
writetable(St,fullfile(OutDir,'FilamentCountPerSetting.csv'));
disp(St)
save(fullfile(OutDir,'FilamentCounter.mat'),'Ccell','VoltageLag','Tp','Fq','UseQuanta','Zt','Kappa','Cd', ...
    'KSigma','KDet','Eta','HVRatio','ProbeGain','IsoSeg','IsoQ','Pre','Post','QRange');

%% Section 9 - Figures
Col = lines(numel(Settings));

% Cell capacitance check
fig = figure('Color','w','Position',[40 40 1500 450]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile
plot(LagScan*1e9,(CalRms./min(CalRms,[],2))','Color',[.7 .7 .7]); hold on
plot(LagScan*1e9,mean(CalRms./min(CalRms,[],2),1),'k','LineWidth',1.5); xline(VoltageLag*1e9,'r--');
xlabel('assumed voltage-channel lag / ns'); ylabel('rms residual / its minimum'); box on
title(sprintf('Probe lag from %d discharge-free captures: %.0f ns',height(Cal),VoltageLag*1e9))
nexttile
errorbar(mean(CalV,1)/1e3,mean(CalC,1)*1e12,std(CalC,0,1)*1e12,'ko-','MarkerFaceColor','k'); hold on
yline(Ccell*1e12,'r--',sprintf('C_{cell} = %.2f pF',Ccell*1e12)); yline(7.799,'b:','old 7.799 pF');
xlabel('peak |V| in the fit window / kV'); ylabel('C / pF'); box on
title('C_{cell} by window (duty 0.08, no discharge)')
nexttile
scatter(Caps.Setting,Caps.CQuiet_pF,25,'filled'); hold on; yline(Ccell*1e12,'r--');
xlabel('duty setting'); ylabel('C on discharge-free stretches / pF'); box on
title({'Same fit inside every capture','(above C_{cell}: hidden discharge current)'})
SaveFigure(fig,OutDir,'CellCapacitanceCheck')

% Pulse template
fig = figure('Color','w','Position',[40 40 1300 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
S = IsoSeg(Use,:); S = S./max(S,[],2);
plot(((0:size(S,2)-1) - Pre)*dt*1e9,S','Color',[.85 .85 .85]); hold on
plot(Tp.tFine*1e9,Tp.h/max(Tp.h),'k','LineWidth',1.5); xlim([-40 400]); ylim([-1 1.2]); box on
xlabel('t - t_{peak} / ns'); ylabel('current / peak')
title(sprintf('Single-pulse response: %d isolated pulses, FWHM %.0f ns',Tp.N,TpFWHM*1e9))
nexttile
plot(Tp.tFine*1e9,cumtrapz(Tp.h)*dt/Tp.Up,'k','LineWidth',1.5); xlim([-40 800]); yline(1,':'); box on
xlabel('t - t_{peak} / ns'); ylabel('cumulative share of the pulse charge')
title('The charge is only complete once the ringing has died out')
SaveFigure(fig,OutDir,'PulseTemplate')

% Example capture: events, pulses and a few decomposed events
f = find(strcmp(Caps.Name,ExampleFile),1);
if ~isempty(f)
    [t,V,I] = ReadScopeCSV(Caps.File{f}); dt = t(2)-t(1);
    V = ShiftSamples(V*HVRatio,-round(VoltageLag/dt)); I = I*ProbeGain;
    Ev = DetectDischargeEvents(t,V,I,Ccell,struct('KSigma',KSigma));
    E = Events(Events.Capture == f,:); P = Pulses(Pulses.Capture == f,:);
    fig = figure('Color','w','Position',[40 40 1500 900]);
    tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
    nexttile([1 3])
    us = (t - Ev.tSwitch)*1e6; w = us > -2 & us < 120;
    plot(us(w),Ev.Ic(w),'Color',[.6 .6 .6]); hold on
    yline([Ev.Threshold -Ev.Threshold],'b:');
    stem(P.t_us,P.q_nC/20,'r','Marker','none');
    plot(us(w),V(w)/1e5,'Color',[0.85 0.5 0.1]);
    xlim([-2 120]); box on; xlabel('t - t_{switch} / \mus'); ylabel('I / A')
    title(sprintf('%s: %d events, %d pulses (red: charge/20, nC per A); orange: V/100 kV', ...
        strrep(ExampleFile,'_','\_'),height(E),height(P)))
    [~,o] = sort(abs(E.Q_nC),'descend');
    for j = o(round(linspace(1,numel(o),min(6,numel(o)))))'
        nexttile
        e = find(Events.Capture == f & Events.Start == E.Start(j),1);
        [Comp,Fit] = DecomposeEvent(Segs{e},Tp,Caps.Sigma_mA(f)/1e3,struct('KDet',KDet,'Eta',Eta));
        tt = (0:numel(Segs{e})-1)*dt*1e9;
        plot(tt,Segs{e},'Color',[.6 .6 .6]); hold on; plot(tt,Fit,'r');
        if ~isempty(Comp), xline(Comp.tau_ns,'b:'); end
        xlim([0 min(tt(end),500)]); box on
        xlabel('ns from event start'); ylabel('pol \times I / A')
        title(sprintf('%.2f nC, %d pulse(s)',abs(E.Q_nC(j)),height(Comp)))
    end
    SaveFigure(fig,OutDir,'EventDecomposition')
end

% Charge quanta
fig = figure('Color','w','Position',[40 40 1300 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
x = abs(Events.Q_nC);
histogram(x,0:0.1:8,'Normalization','count','FaceColor',[.7 .7 .7]); hold on
xx = linspace(QRange(1),QRange(2),800)'; Bw = 0.1*Fq.NFit;
plot(xx,Bw*QuantaDensity(xx,Fq),'r','LineWidth',1.3); plot(xx,Bw*SmoothDensity(xx,Fq),'b','LineWidth',1.3);
xline(QRange,'k--','HandleVisibility','off'); xlim([0 8]); box on
xlabel('event charge |Q| / nC'); ylabel('events per 0.1 nC')
legend('events',sprintf('best quantized fit, q_1 = %.2f nC',Fq.q1),'smooth (no humps)')
title(sprintf('All events, all settings (N = %d; fit range dashed)',height(Events)))
nexttile
histogram(Fq.GainBoot,20,'FaceColor',[.7 .7 .7]); hold on
xline(Fq.LogLQuanta - Fq.LogLSmooth,'r','LineWidth',2); box on
xlabel('log-likelihood gain of the quantized model'); ylabel('bootstrap replicates')
title(sprintf('Smooth-model replicates vs the data (red): p = %.2f',Fq.PBoot))
SaveFigure(fig,OutDir,'ChargeQuanta')

% Lissajous: shot-averaged, drift-corrected measured loop vs reconstruction from events,
% and per half-cycle the Lissajous charge against the event charge
fig = figure('Color','w','Position',[40 40 1700 800]);
tiledlayout(2,numel(Settings),'TileSpacing','compact','Padding','compact');
tg = (-5:0.01:150)*1e-6;
for s = 1:numel(Settings)
    nexttile(s); hold on
    c = find(Caps.Setting == SetVal(s));
    [Vm,Qm,Qr] = deal(zeros(numel(c),numel(tg)));
    for i = 1:numel(c)
        W = Wave{c(i)};
        Vm(i,:) = interp1(W.t,W.V,tg);
        Qm(i,:) = Ccell*Vm(i,:) + interp1(W.t,Rc{c(i)},tg);
        E = Events(Events.Capture == c(i),:);
        Qr(i,:) = Ccell*Vm(i,:) + arrayfun(@(x) sum(E.Q_nC(E.tPeak_us*1e-6 < x)),tg)*1e-9;
    end
    plot(mean(Vm,1)/1e3,mean(Qm,1)*1e9,'Color',[.5 .5 .5]);
    plot(mean(Vm,1)/1e3,mean(Qr,1)*1e9,'r');
    box on; xlabel('V / kV'); ylabel('Q / nC')
    title(sprintf('duty %s, mean of %d shots',Settings{s},numel(c)))
    if s == 1, legend('measured \int I dt, drift-corrected','C_{cell}V + event charges','Location','northwest'); end
    nexttile(numel(Settings) + s); hold on
    Hs = HC(HC.Setting == SetVal(s) & ~isnan(HC.QMeas_nC),:);
    Hs.K = round(Hs.tPeak_us/(median(diff(HC.tPeak_us(HC.Capture == c(1))))));
    G = findgroups(Hs.K);
    Vk = splitapply(@(x) mean(abs(x)),Hs.Vpk_kV,G);
    Mm = splitapply(@mean,Hs.Pol.*Hs.QMeas_nC,G); Ms = splitapply(@(x) std(x)/sqrt(numel(x)),Hs.Pol.*Hs.QMeas_nC,G);
    Pm = splitapply(@mean,Hs.Pol.*Hs.QApp_nC,G);
    errorbar(Vk,Mm,Ms,'ko','MarkerFaceColor',[.5 .5 .5]); plot(Vk,Pm,'rs','MarkerFaceColor','r');
    yline(0,':'); box on; xlabel('half-cycle peak |V| / kV'); ylabel('discharge charge per half-cycle / nC')
    if s == 1, legend('Lissajous (mean \pm s.e.)','detected events','Location','northwest'); end
end
SaveFigure(fig,OutDir,'Lissajous')

% Voltage sweep: per half-cycle against its peak voltage
fig = figure('Color','w','Position',[40 40 1500 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
Mk = {'v','^'};
Panels = {'Filaments','QApp_nC','QPerFil_nC','VOnset_kV','UgOnset_kV','Energy_uJ'};
YLab = {'filaments per half-cycle','event charge per half-cycle / nC','charge per filament / nC', ...
    'applied voltage at onset / kV','gap voltage at onset / kV','energy per half-cycle / \muJ'};
for p = 1:numel(Panels)
    nexttile; hold on
    for s = 1:numel(Settings)
        for pp = 1:2
            r = HC.Setting == SetVal(s) & HC.Pol == 2*pp - 3;
            y = HC.(Panels{p})(r);
            if ismember(Panels{p},{'QApp_nC','VOnset_kV','UgOnset_kV'}), y = HC.Pol(r).*y; end   % own polarity > 0
            h = scatter(abs(HC.Vpk_kV(r)),y,14,Col(s,:),Mk{pp},'filled','MarkerFaceAlpha',0.6);
            if pp == 1, h.HandleVisibility = 'off'; end
        end
    end
    box on; xlabel('half-cycle peak |V| / kV'); ylabel(YLab{p})
    if p == 1, legend("duty " + string(Settings),'Location','northwest'); end
    if ~GapOK && ismember(Panels{p},{'UgOnset_kV','Energy_uJ'})
        text(0.5,0.5,{'needs the dielectric','(EpsR, DielThick\_mm)'},'Units','normalized','HorizontalAlignment','center')
    end
end
sgtitle('Ring-down as a voltage sweep (\nabla negative, \Delta positive half-cycles; voltages and charges signed so their own polarity is > 0)')
SaveFigure(fig,OutDir,'VoltageSweep')

% Counts and charge per filament by setting
fig = figure('Color','w','Position',[40 40 1400 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
b = bar(1:height(St),[St.PosFilaments St.NegFilaments St.Filaments St.PulsesExploratory]);
hold on; errorbar(b(3).XEndPoints,St.Filaments,St.FilamentsSD,'k.','HandleVisibility','off')
set(gca,'XTick',1:height(St),'XTickLabel',St.Setting); box on
ylabel('per capture (mean of shots)'); xlabel('duty setting')
legend('positive','negative','all filaments (\pm s.d.)','decomposed pulses (exploratory)','Location','northwest')
title({'Filaments per excitation (lower bound)',CountRule})
nexttile; hold on
for s = 1:numel(Settings)
    E = Events(Events.Setting == SetVal(s),:);
    if isempty(E), continue, end
    histogram(abs(E.Q_nC)./E.N,10.^(-1:0.1:1.3),'DisplayStyle','stairs','EdgeColor',Col(s,:),'LineWidth',1.3, ...
        'DisplayName',"duty " + Settings{s});
end
box on; set(gca,'XScale','log'); xlim([0.1 16]); xlabel('apparent charge per filament / nC'); ylabel('filaments'); legend
title('Charge per filament (detected events only)')
SaveFigure(fig,OutDir,'FilamentCounts')

%% Local functions
function T = ListCaptures(Dir,Pattern,Settings)
% One row per capture of the given settings, in shot order.
T = table();
for s = 1:numel(Settings)
    d = dir(fullfile(Dir,sprintf(Pattern,Settings{s})));
    [~,o] = sort({d.name}); d = d(o);
    for i = 1:numel(d)
        tok = regexp(d(i).name,'_([\d.]+)_(\d+)_Analog','tokens','once');
        T = [T; table({fullfile(d(i).folder,d(i).name)},{d(i).name},str2double(tok{1}),str2double(tok{2}), ...
            'VariableNames',{'File','Name','Setting','Shot'})]; %#ok<AGROW>
    end
end
end

function d = SmoothDeriv(V,dt)
% dV/dt with the same 1.6 us Savitzky-Golay smoothing as DetectDischargeEvents.
d = smooth(gradient(V)./dt,2*floor(round(1.6e-6/dt)/2)+1,'sgolay',3)';
end

function i = SwitchIndex(t,V)
% First sample where V leaves its starting level by more than 300 V.
i = find(abs(V - median(V(t < t(1) + 20e-6))) > 300,1);
end

function y = ShiftSamples(x,k)
% y(i) = x(i - k): positive k delays the trace; the ends are held.
n = numel(x); i = min(max((1:n) - k,1),n);
y = x(i);
end

function [H,Rc] = SplitHalfCycles(W,Ccell,tEv)
% Split a ring-down at its voltage peaks and measure each half-cycle's discharge charge from
% the integrated current. After a voltage peak no discharge can start until |V| has fallen
% well back, so on that stretch (peak + 0.5 us until |V| < 30% of the peak, cut short at any
% event) Q - Ccell*V is flat apart from the CH2 baseline drift. A spline drift (knots every
% 40 us) plus one level per stretch is fitted to those stretches; the step between
% consecutive levels is the discharge charge of the half-cycle in between.
%   W: t (s after the switch), V, Q; tEv: [start stop] times of the capture's events (s).
% Returns H (table: tPeak_us, Vpk_kV, Pol, tWin_us = [start stop] of the discharge window,
% QMeas_nC) and Rc, the drift-corrected Q - Ccell*V on W.t.
dtw = W.t(2) - W.t(1);
Vs = movmean(W.V,max(1,round(0.5e-6/dtw)));
[~,ip] = findpeaks(abs(Vs),'MinPeakDistance',round(12e-6/dtw),'MinPeakHeight',1000);
ip = ip(W.t(ip) > 0 & W.t(ip) < W.t(end) - 10e-6);
r = W.Q - Ccell*W.V;
r = r - mean(r(W.t > -4.5e-6 & W.t < -2.5e-6));
Seg = zeros(numel(ip)+1,2);
Seg(1,:) = [find(W.t > -4.5e-6,1) find(W.t < -2.5e-6,1,'last')];
for k = 1:numel(ip)
    a = ip(k) + round(0.5e-6/dtw);
    b = a - 1 + find(abs(Vs(a:end)) < 0.3*abs(Vs(ip(k))),1);
    if isempty(b), b = numel(W.t); end
    ev = tEv(tEv(:,2) >= W.t(a) & tEv(:,1) <= W.t(b),:);
    if ~isempty(ev), b = find(W.t < min(ev(:,1)) - 1e-6,1,'last'); end
    Seg(k+1,:) = [a b];
end
ok = Seg(:,2) - Seg(:,1) >= round(1e-6/dtw);
idx = []; lab = [];
for k = find(ok)'
    s = Seg(k,1):5:Seg(k,2); idx = [idx s]; lab = [lab k*ones(size(s))]; %#ok<AGROW>
end
T0 = W.t(Seg(1,1)); T1 = W.t(max(Seg(ok,2)));
Knots = linspace(T0,T1,max(2,round((T1 - T0)/40e-6)) + 1);
B = SplineBasis(Knots,W.t(idx));
Lm = full(sparse(1:numel(idx),lab,1,numel(idx),size(Seg,1)));
x = [Lm(:,ok) B(:,2:end)] \ r(idx)';
L = nan(size(Seg,1),1); L(ok) = x(1:nnz(ok));
D = SplineBasis(Knots,min(max(W.t,T0),T1))*[0; x(nnz(ok)+1:end)];
Rc = r - D';
n = numel(ip);
tWin = [W.t(Seg(1:n,2))' W.t(Seg(2:n+1,1))']*1e6;
H = table(W.t(ip)'*1e6,W.V(ip)'/1e3,sign(W.V(ip))',tWin,(L(2:n+1) - L(1:n))*1e9, ...
    'VariableNames',{'tPeak_us','Vpk_kV','Pol','tWin_us','QMeas_nC'});
end

function B = SplineBasis(Knots,tq)
% Columns: the not-a-knot cubic spline through a unit value at one knot, at times tq.
B = zeros(numel(tq),numel(Knots));
for q = 1:numel(Knots)
    y = zeros(size(Knots)); y(q) = 1;
    B(:,q) = spline(Knots,y,tq(:));
end
end

function f = QuantaDensity(x,Fq)
nn = 1:Fq.NMax; P = Fq.ParQuanta;
w = P(4).^nn./factorial(nn); w = w/sum(w);
mu = nn*P(1); sd = sqrt(P(2)^2 + nn*P(3)^2);
R = Fq.QRange;
f = sum(w.*normpdf(x,mu,sd),2)/sum(w.*(normcdf(R(2),mu,sd) - normcdf(R(1),mu,sd)));
end

function f = SmoothDensity(x,Fq)
z = Fq.ParSmooth; w = 1/(1 + exp(-z(5))); R = Fq.QRange;
s1 = exp(z(2)); s2 = exp(z(4));
f = w*lognpdf(x,z(1),s1) + (1 - w)*lognpdf(x,z(3),s2);
f = f/(w*(logncdf(R(2),z(1),s1) - logncdf(R(1),z(1),s1)) + (1 - w)*(logncdf(R(2),z(3),s2) - logncdf(R(1),z(3),s2)));
end

function SaveFigure(fig,Folder,Name)
% Save a figure as .png and .fig and leave it open, titled with its file name.
fig.Name = Name;
fig.Visible = 'on';
exportgraphics(fig,fullfile(Folder,[Name '.png']))
savefig(fig,fullfile(Folder,[Name '.fig']))
end
