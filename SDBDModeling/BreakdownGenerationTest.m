%% Breakdown-generation test on the charge derivative
%
% Hypothesis: the peaks in the positive half-wave dQp/dUr (paper Fig 6) are
% successive breakdown "generations". A positive surface discharge charges the
% dielectric until the field at the triple junction is cancelled; the sites go
% quiet until Ur has risen by roughly the burning voltage Ug again, then fire
% together. If that is right, across applied amplitudes:
%   1. the peaks sit at a fixed voltage swing since the reversal (the previous
%      Ur extremum), dU = |Ur - Ur_extremum|, and so drift in Ur with slope -1
%      against the amplitude. A fixed instantaneous-voltage threshold would
%      instead give slope 0 in Ur. (In the classic Manley picture ignition
%      itself comes a fixed swing of about 2 Ug after the reversal.)
%   2. consecutive peaks are spaced by a fixed swing, of the order of Ug;
%   3. peak heights stay about constant (paper: ~200 and ~420 pF).
%
% For every trace of a voltage sweep: Resr and Coff from the plasma-off traces,
% Ur = V - I*Resr, phase-averaged Q and Ur, fourth-order Savitzky-Golay
% derivative, turning points masked (paper Section II-A). dQp/dUr is the total
% dQ/dUr minus its off-phase level on each branch, so it reads 0 before
% ignition. Ignition = first bin where dQp/dUr reaches IgnFrac of the branch
% maximum and stays there. Ug is estimated per trace as half the ignition swing
% (Manley: the gap goes from -Ug to +Ug while the dielectric charge holds).
% Peaks of the high-amplitude traces are grouped into bands by k-means on their
% swing since reversal, and every peak is assigned to the nearest band.
%
% Section 5 then assigns the generations one trace at a time, gated in
% dQp/dUr: generation g is the first peak, at least MinGenGap after the
% previous accepted event, whose height lies inside GenWindows(g,:). A stray or
% missed peak therefore cannot renumber the ones after it. Onset, generation 1,
% 2 and 3 are plotted as the swing since the previous half-wave's extremum,
% which should give a staircase.
%
% Resr: the off-trace Re(Z) includes the dielectric loss tan(d)/(w*Coff), which
% belongs to the device, not to the series resistance subtracted from V. With
% the Hioki tan(d) of the device (SDBDCleanVersion.m) it is removed:
% Resr = Re(Z) - tan(d)/(w*Coff).
%
% Outputs: Output/GenerationTest (figures left open, .png and .fig).

clear

%% Section 1 - Options
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
VC       = fullfile(RepoRoot,'Data','SDBDVoltageandCurrent','VoltageandCurrent','2025');
Sweeps(1).Name = '1600 Hz';
Sweeps(1).On   = fullfile(VC,'1600 Hz','bighex','bighex','vsweep','PlasmaOn');
Sweeps(1).Off  = fullfile(VC,'1600 Hz','bighex','bighex','vsweep','PlasmaOff');
Sweeps(1).F0   = 1600;
Sweeps(2).Name = '3000 Hz';
Sweeps(2).On   = fullfile(VC,'3000 Hz','bighexvsweep','bighexvsweep','PlasmaOn');
Sweeps(2).Off  = fullfile(VC,'3000 Hz','bighexvsweep','bighexvsweep','PlasmaOff');
Sweeps(2).F0   = 3000;
OutDir = fullfile(RepoRoot,'Output','GenerationTest');
if ~exist(OutDir,'dir'), mkdir(OutDir); end

AmpGain          = 500;      % V/V, Trek voltage monitor
CurrentGain      = 6e-3;     % A/V, Trek current monitor
PhaseBins        = 1000;
SGOrder          = 4;        % paper: fourth-degree Savitzky-Golay derivative
SGWindow         = 41;
TurningPointFrac = 0.10;     % paper: mask where dUr/dphase < 10% of its maximum
BaselineFrac     = 0.15;     % first share of each branch (after the mask) taken as off-phase
IgnFrac          = 0.30;     % ignition = first bin at this share of the branch maximum...
IgnHoldBins      = 15;       % ...that stays above half of it for this many bins
NBands           = 3;        % peak bands to track on the positive branch
BandMinAmp       = 2000;     % V; traces used to define the bands
MinBranchPeak    = 60;       % pF; traces whose branch never gets this high are below breakdown
PeakProminence   = 25;       % pF, for the peaks in dQp/dUr
ResrOverride     = [];       % [R1600 R3000] in Ohm to replace the off-trace Resr; NaN keeps it
LossTangent      = 0.019;    % device tan(d) from the Hioki (SDBDCleanVersion.m); 0 = paper's Re(Z)
IgnOnsetFrac     = 0.10;     % discharge onset = first bin at this share of the branch maximum (held)
ImagingCache     = fullfile(RepoRoot,'Output','Cache','BigHexDerived.mat');   % ICCD overlay (3 kHz imaging day)
OnsetMinPeak     = 150;      % pF; onset only reported when the branch reaches this (below, noise sets it)
MinGenGap        = 500;      % V of swing between accepted generations
Sweeps(1).GenWindows = [140 260; 300 480; 430 700];   % pF, dQp/dUr window of generations 1-3
Sweeps(2).GenWindows = [150 290; 280 460; 380 650];

%% Section 2 - Every trace of every sweep
Res = table;
Curves = struct([]);
for s = 1:numel(Sweeps)
    Sw = Sweeps(s);
    [ReZ,Coff] = OffMeasurement(Sw.Off,AmpGain,CurrentGain,Sw.F0);
    Resr = ReZ - LossTangent/(2*pi*Sw.F0*Coff);         % external series resistance only
    if ~isempty(ResrOverride) && ~isnan(ResrOverride(s)), Resr = ResrOverride(s); end
    fprintf('\n%s: Re(Z) = %.2f kOhm, dielectric loss %.2f kOhm -> Resr = %.2f kOhm; Coff = %.1f pF\n', ...
            Sw.Name,ReZ/1e3,LossTangent/(2*pi*Sw.F0*Coff)/1e3,Resr/1e3,Coff*1e12);
    Files = dir(fullfile(Sw.On,'*.csv'));
    [~,ix] = sort({Files.name}); Files = Files(ix);
    for n = 1:numel(Files)
        [t,V,I] = ReadScope(fullfile(Files(n).folder,Files(n).name),AmpGain,CurrentGain);
        f = FitFrequency(t,V,Sw.F0);
        P = 1/f;
        Ncyc = floor((t(end) - t(1))/P);
        if Ncyc < 2, continue, end
        k = t < t(1) + Ncyc*P;
        t = t(k);  V = V(k);  I = I(k) - mean(I(k));
        Ur = V - I*Resr;
        Q  = cumtrapz(t,I);  Q = Q - mean(Q);
        b   = min(floor(mod((t - t(1))/P,1)*PhaseBins) + 1,PhaseBins);
        cnt = accumarray(b,1,[PhaseBins 1]);
        Ub  = accumarray(b,Ur,[PhaseBins 1])./cnt;
        Qb  = accumarray(b,Q, [PhaseBins 1])./cnt;
        [~,g] = sgolay(SGOrder,SGWindow);
        dU = CircDeriv(Ub,g);
        Mask = abs(dU) < TurningPointFrac*max(abs(dU));
        Ctot = CircDeriv(Qb,g)./dU;
        Ctot(Mask) = NaN;
        Amp = max(abs(Ub));
        [~,iMin] = min(Ub);  [~,iMax] = max(Ub);
        % Branches in order: positive = Ur rising (min -> max), negative = falling
        Br{1} = mod((iMin-1):(iMin-1+mod(iMax-iMin,PhaseBins)),PhaseBins)' + 1;
        Br{2} = mod((iMax-1):(iMax-1+mod(iMin-iMax,PhaseBins)),PhaseBins)' + 1;
        for p = 1:2
            idx = Br{p};
            idx = idx(~isnan(Ctot(idx)));
            nb  = max(3,round(BaselineFrac*numel(idx)));
            Base = median(Ctot(idx(1:nb)));
            C  = (Ctot(idx) - Base)*1e12;                   % dQp/dUr, pF
            U  = Ub(idx);
            Urev = Ub(Br{p}(1));                            % the previous extremum
            dUs  = abs(U - Urev);                           % swing since reversal
            if max(C) < MinBranchPeak, continue, end
            % Ignition: first bin at IgnFrac of the maximum that stays above half of it
            Above = C >= IgnFrac*max(C);
            Hold  = movmin(double(C >= 0.5*IgnFrac*max(C)),[0 IgnHoldBins-1]) > 0;
            j0 = find(Above & Hold,1,'first');
            if isempty(j0), continue, end
            Uign = U(j0);  SwIgn = dUs(j0);
            % Onset: where dQp/dUr first lifts off zero and stays up (start of generation 1)
            HoldO = movmin(double(C >= 0.5*IgnOnsetFrac*max(C)),[0 IgnHoldBins-1]) > 0;
            jO = find(C >= IgnOnsetFrac*max(C) & HoldO,1,'first');
            if isempty(jO), jO = j0; end
            SwOnset = dUs(jO);
            if max(C) < OnsetMinPeak, SwOnset = NaN; end
            Ug   = SwIgn/2;                                 % Manley: -Ug -> +Ug at fixed dielectric charge
            On   = (1:numel(C))' >= j0;
            [pk,loc] = findpeaks(C(On),'MinPeakProminence',PeakProminence);
            loc = loc + j0 - 1;
            for q = 1:numel(pk)
                Res = [Res; table(string(Sw.Name),string(Files(n).name),p,Amp,Uign,SwIgn,Ug,q, ...
                       pk(q),U(loc(q)),dUs(loc(q)),'VariableNames',{'Sweep','File','Branch', ...
                       'UrAmplitude_V','UrIgnition_V','IgnitionSwing_V','Ug_V','PeakNo','Peak_pF', ...
                       'PeakUr_V','PeakSwing_V'})]; %#ok<AGROW>
            end
            Curves(end+1).Sweep = s; %#ok<SAGROW>
            Curves(end).Branch = p;  Curves(end).Amp = Amp;  Curves(end).U = U;  Curves(end).C = C;
            Curves(end).dU = dUs;    Curves(end).Uign = Uign; Curves(end).Ug = Ug; Curves(end).SwIgn = SwIgn;
            Curves(end).SwOnset = SwOnset;  Curves(end).File = string(Files(n).name);
        end
    end
end
writetable(Res,fullfile(OutDir,'GenerationPeaks.csv'));

%% Section 3 - Are the peaks locked to the reversal, or to the instantaneous Ur?
% For each band: position against Ur amplitude, as a straight line. Locked to
% the reversal -> slope -1 in Ur and 0 in swing; a fixed Ur threshold -> slope 0
% in Ur and +1 in swing.
BrName = {'positive','negative'};
Res.Band = nan(height(Res),1);
Bands = cell(numel(Sweeps),2);
rng(3)
for s = 1:numel(Sweeps)
    for p = 1:2
        sel = Res.Sweep == Sweeps(s).Name & Res.Branch == p;
        c = Curves([Curves.Sweep] == s & [Curves.Branch] == p);
        if isempty(c), continue, end
        fprintf('\n%s, %s half-wave: %d traces above breakdown, Ur amplitude %.0f-%.0f V\n', ...
                Sweeps(s).Name,BrName{p},numel(c),min([c.Amp]),max([c.Amp]));
        Hi = [c.Amp] >= BandMinAmp;
        fprintf('  ignition swing since reversal %.0f +- %.0f V (amplitude >= %d V), so Ug ~ %.0f V; peaks per trace %.1f\n', ...
                mean([c(Hi).SwIgn]),std([c(Hi).SwIgn]),BandMinAmp,mean([c(Hi).Ug]),nnz(sel)/numel(c));
        [bI,ciI] = SlopeCI([c.Amp]',[c.SwIgn]');
        fprintf('  ignition swing vs amplitude: slope %+.2f [%+.2f, %+.2f] V/V\n',bI,ciI(1),ciI(2));
        r = Res(sel & Res.UrAmplitude_V >= BandMinAmp,:);
        nb = NBands;
        if p == 2, nb = min(NBands,2); end
        if height(r) < 3*nb, continue, end
        [~,Cent] = kmeans(r.PeakSwing_V,nb,'Replicates',20);
        Cent = sort(Cent);
        Bands{s,p} = Cent;
        [~,Res.Band(sel)] = min(abs(Res.PeakSwing_V(sel) - Cent'),[],2);
        for q = 1:nb
            rq = Res(sel & Res.Band == q & Res.UrAmplitude_V >= 1500,:);
            if height(rq) < 5, continue, end
            [bU,ciU] = SlopeCI(rq.UrAmplitude_V,rq.PeakUr_V);
            [bS,ciS] = SlopeCI(rq.UrAmplitude_V,rq.PeakSwing_V);
            fprintf(['  band %d (%d peaks): swing %.0f +- %.0f V, Ur %.0f +- %.0f V, height %.0f +- %.0f pF\n' ...
                     '      slope vs amplitude: in Ur %+.2f [%+.2f, %+.2f], in swing %+.2f [%+.2f, %+.2f]\n'], ...
                    q,height(rq),mean(rq.PeakSwing_V),std(rq.PeakSwing_V),mean(rq.PeakUr_V),std(rq.PeakUr_V), ...
                    mean(rq.Peak_pF),std(rq.Peak_pF),bU,ciU(1),ciU(2),bS,ciS(1),ciS(2));
        end
        if numel(Cent) >= 2
            fprintf('  band spacing (swing): %s V; compare Ug ~ %.0f V\n', ...
                    sprintf('%.0f ',diff(Cent)),mean([c(Hi).Ug]));
        end
    end
end
writetable(Res,fullfile(OutDir,'GenerationPeaks.csv'));

%% Section 4 - Figures, one per sweep
for s = 1:numel(Sweeps)
    fig = figure('Color','w','Position',[40 40 1500 900]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    c = Curves([Curves.Sweep] == s);
    if isempty(c), continue, end
    AmpAll = [c.Amp];  Cm = parula(256);
    ColOf = @(a) Cm(1 + round(255*(a - min(AmpAll))/max(eps,max(AmpAll) - min(AmpAll))),:);

    nexttile; hold on   % Fig 6 as published: against normalised Ur
    for i = 1:numel(c), plot(c(i).U/c(i).Amp,c(i).C,'Color',ColOf(c(i).Amp)); end
    xlabel('U_r / max(U_r)'); ylabel('dQ_p/dU_r / pF'); title([Sweeps(s).Name ': as in Fig 6']); box on; grid on
    colormap(gca,Cm); clim([min(AmpAll) max(AmpAll)]); cb = colorbar; cb.Label.String = 'U_r amplitude / V';

    for p = 1:2
        nexttile; hold on
        cp = c([c.Branch] == p);
        for i = 1:numel(cp), plot(cp(i).dU,cp(i).C,'Color',ColOf(cp(i).Amp)); end
        r = Res(Res.Sweep == Sweeps(s).Name & Res.Branch == p,:);
        scatter(r.PeakSwing_V,r.Peak_pF,14,'k','filled','MarkerFaceAlpha',0.5)
        if ~isempty(Bands{s,p}), xline(Bands{s,p},'r--'); end
        xlabel('Swing since reversal |U_r - U_{r,extremum}| / V'); ylabel('dQ_p/dU_r / pF'); box on; grid on
        title(sprintf('%s half-wave, aligned at the reversal',[upper(BrName{p}(1)) BrName{p}(2:end)]))
        colormap(gca,Cm); clim([min(AmpAll) max(AmpAll)])
    end

    r = Res(Res.Sweep == Sweeps(s).Name & Res.Branch == 1,:);
    Mk = {'o','s','^'};
    nexttile; hold on
    for q = 1:3
        rq = r(r.Band == q,:);
        scatter(rq.UrAmplitude_V,rq.PeakUr_V,24,Mk{q},'filled','MarkerFaceAlpha',0.6)
    end
    xlabel('U_r amplitude / V'); ylabel('Peak position U_r / V'); box on; grid on
    legend('Band 1','Band 2','Band 3','Location','best'); title('Positive peaks in U_r (locked to reversal: slope -1)')
    nexttile; hold on
    for q = 1:3
        rq = r(r.Band == q,:);
        scatter(rq.UrAmplitude_V,rq.PeakSwing_V,24,Mk{q},'filled','MarkerFaceAlpha',0.6)
    end
    cp = c([c.Branch] == 1);
    plot([cp.Amp],[cp.SwIgn],'k.','MarkerSize',10)
    xlabel('U_r amplitude / V'); ylabel('Swing since reversal / V'); box on; grid on
    legend('Band 1','Band 2','Band 3','Ignition','Location','best'); title('Positive peaks in swing since reversal (locked: flat)')
    nexttile; hold on
    for q = 1:3
        rq = r(r.Band == q,:);
        scatter(rq.UrAmplitude_V,rq.Peak_pF,24,Mk{q},'filled','MarkerFaceAlpha',0.6)
    end
    xlabel('U_r amplitude / V'); ylabel('Peak height / pF'); box on; grid on
    legend('Band 1','Band 2','Band 3','Location','best'); title('Positive peak heights')
    SaveFigure(fig,OutDir,['GenerationTest_' strrep(Sweeps(s).Name,' ','')])
end

%% Section 5 - Generations gated in dQp/dUr, as a staircase in swing since reversal
Gen = table;
for s = 1:numel(Sweeps)
    W  = Sweeps(s).GenWindows;
    c  = Curves([Curves.Sweep] == s & [Curves.Branch] == 1);
    rS = Res(Res.Sweep == Sweeps(s).Name & Res.Branch == 1,:);
    for i = 1:numel(c)
        rp = sortrows(rS(rS.File == c(i).File,:),'PeakSwing_V');
        Sw = nan(1,size(W,1));  Ht = Sw;
        Prev = c(i).SwOnset;  Gap = 0;  used = false(height(rp),1);
        for g = 1:size(W,1)
            k = find(~used & rp.PeakSwing_V >= Prev + Gap & rp.Peak_pF >= W(g,1) & rp.Peak_pF <= W(g,2),1,'first');
            if isempty(k), continue, end
            Sw(g) = rp.PeakSwing_V(k);  Ht(g) = rp.Peak_pF(k);
            used(1:k) = true;  Prev = Sw(g);  Gap = MinGenGap;
        end
        Gen = [Gen; table(string(Sweeps(s).Name),c(i).File,c(i).Amp,c(i).SwOnset,Sw(1),Sw(2),Sw(3),Ht(1),Ht(2),Ht(3), ...
               'VariableNames',{'Sweep','File','UrAmplitude_V','Onset_V','Gen1_V','Gen2_V','Gen3_V', ...
               'Gen1_pF','Gen2_pF','Gen3_pF'})]; %#ok<AGROW>
    end
end
writetable(Gen,fullfile(OutDir,'GenerationStaircase.csv'));

% The staircase starts at the reversal (swing 0). The onset is kept in the CSV
% and printout only: it marks particularly early ignition regions, not a new
% generation.
EvName = {'Reversal','Generation 1','Generation 2','Generation 3'};

% ICCD overlay: positive half-wave gates of the 3 kHz imaging day, as swing
% since the previous Ur extremum, and the peaks of the image intensity
Img = struct('Ok',false);
if exist(ImagingCache,'file')
    Dc = load(ImagingCache,'UrGate','HalfWave','HumpStart','HumpEnd','Tot337','SideProf','SideLit','DepthMM');
    h  = find(Dc.HalfWave(Dc.HumpStart) > 0,1);
    gI = Dc.HumpStart(h):Dc.HumpEnd(h);
    Img.Ok    = true;
    Img.Amp   = -min(Dc.UrGate);                    % previous extremum (Ur minimum)
    Img.Swing = Dc.UrGate(gI) - min(Dc.UrGate);
    Img.Tot   = Dc.Tot337(gI)/max(Dc.Tot337(gI));
    near  = Dc.DepthMM >= 0 & Dc.DepthMM <= 0.25;
    Edge  = squeeze(mean(Dc.SideProf(:,near,gI,2),2));
    Lit   = Dc.SideLit(:,gI) >= 0.5;
    WasLit = [false(size(Lit,1),1) Lit(:,1:end-1)];
    Img.Edge = arrayfun(@(n) mean(Edge(WasLit(:,n),n)),1:numel(gI));
    Img.Edge = Img.Edge/max(Img.Edge);
    [~,l1] = findpeaks(Img.Tot,'MinPeakProminence',0.03);
    [~,l2] = findpeaks(Img.Edge,'MinPeakProminence',0.03);
    Img.PkTot  = Img.Swing(l1);
    Img.PkEdge = Img.Swing(l2);
    fprintf('\nICCD (3 kHz imaging day, Ur amplitude %.0f V): intensity peaks at swing %s V; edge-brightness peaks at %s V\n', ...
            Img.Amp,sprintf('%.0f ',Img.PkTot),sprintf('%.0f ',Img.PkEdge));
end
EvCol  = [0 0 0; 0.15 0.35 0.85; 0.9 0.45 0.1; 0.85 0.7 0.1];
for s = 1:numel(Sweeps)
    G = Gen(Gen.Sweep == Sweeps(s).Name,:);
    if isempty(G), continue, end
    Ev = [zeros(height(G),1) G.Gen1_V G.Gen2_V G.Gen3_V];
    fprintf('\n%s, positive half-wave staircase (swing since the previous extremum):\n',Sweeps(s).Name);
    okO = ~isnan(G.Onset_V);
    fprintf('  (onset of early ignition regions: median %.0f V, IQR %.0f V, n = %d)\n', ...
            median(G.Onset_V(okO)),iqr(G.Onset_V(okO)),nnz(okO));
    Med = nan(1,4);
    for e = 1:4
        v = Ev(:,e);  ok = ~isnan(v);
        Med(e) = median(v(ok));
        if e > 1 && nnz(ok) >= 5
            b = SlopeCI(G.UrAmplitude_V(ok),v(ok));
            fprintf('  %-13s n = %2d  median %5.0f V  IQR %4.0f V  slope vs amplitude %+.2f', ...
                    EvName{e},nnz(ok),Med(e),iqr(v(ok)),b);
            if e > 1, fprintf('  height %.0f +- %.0f pF',mean(G{ok,['Gen' num2str(e-1) '_pF']}),std(G{ok,['Gen' num2str(e-1) '_pF']})); end
            fprintf('\n');
        end
    end
    fprintf('  steps between medians: %s V\n',sprintf('%.0f ',diff(Med)));
    % Step size from all traces that have a generation, swing = a + step*g
    gg = [];  yy = [];
    for e = 2:4
        ok = ~isnan(Ev(:,e));
        gg = [gg; (e-1)*ones(nnz(ok),1)];  yy = [yy; Ev(ok,e)]; %#ok<AGROW>
    end
    [Step,ciStep] = SlopeCI(gg,yy);
    fprintf('  staircase fit over generations 1-3: step %.0f V [%.0f, %.0f]\n',Step,ciStep(1),ciStep(2));
    FirstStep = Med(2);
    fprintf('  first step (reversal -> generation 1) %.0f V = %.2f x the later step\n',FirstStep,FirstStep/Step);

    fig = figure('Color','w','Position',[40 40 1600 520]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile; hold on
    for e = 1:4
        ok = ~isnan(Ev(:,e));
        scatter(G.UrAmplitude_V(ok),Ev(ok,e),26,EvCol(e,:),'filled','MarkerFaceAlpha',0.7)
    end
    for e = 1:4, yline(Med(e),':','Color',EvCol(e,:),'HandleVisibility','off'); end
    Leg = EvName;
    if Img.Ok
        scatter(Img.Amp*ones(size(Img.PkTot)),Img.PkTot,80,'k^','filled')
        Leg{end+1} = 'ICCD intensity peaks (3 kHz)';
    end
    xlabel('U_r amplitude of the previous half-wave / V'); ylabel('Swing since the previous extremum / V')
    legend(Leg,'Location','northwest'); box on; grid on
    title(sprintf('%s positive half-wave: referenced to the previous extremum',Sweeps(s).Name))

    nexttile; hold on
    for e = 1:4
        ok = ~isnan(Ev(:,e));
        x = (e-1) + 0.12*randn(nnz(ok),1);
        scatter(x,Ev(ok,e),14,EvCol(e,:),'filled','MarkerFaceAlpha',0.35)
    end
    stairs([-0.5 0.5 1.5 2.5 3.5],[Med Med(end)],'k-','LineWidth',2)
    for e = 2:4
        text(e-1.5,(Med(e-1)+Med(e))/2,sprintf('%+.0f V',Med(e)-Med(e-1)),'HorizontalAlignment','center', ...
             'FontWeight','bold','BackgroundColor','w')
    end
    if Img.Ok
        % ICCD intensity peaks, beside the generation each one is nearest to
        for k = 1:numel(Img.PkTot)
            [~,gk] = min(abs(Med(2:4) - Img.PkTot(k)));
            plot(gk + 0.3,Img.PkTot(k),'k^','MarkerFaceColor','k','MarkerSize',9)
            text(gk + 0.42,Img.PkTot(k),sprintf('%.0f V',Img.PkTot(k)),'FontSize',8)
        end
    end
    xticks(0:3); xticklabels(EvName); xlim([-0.6 3.6]); box on; grid on
    ylabel('Swing since the previous extremum / V')
    title(sprintf('Staircase: first step %.0f V (%.2f x), then %.0f V [%.0f, %.0f]',FirstStep,FirstStep/Step,Step,ciStep(1),ciStep(2)))

    nexttile; hold on
    c = Curves([Curves.Sweep] == s & [Curves.Branch] == 1);
    W = Sweeps(s).GenWindows;
    xl = [0 max(arrayfun(@(k) max(c(k).dU),1:numel(c)))];
    for g = 1:size(W,1)
        patch([xl(1) xl(2) xl(2) xl(1)],[W(g,1) W(g,1) W(g,2) W(g,2)],EvCol(g+1,:),'FaceAlpha',0.12, ...
              'EdgeColor','none','HandleVisibility','off')
    end
    for k = 1:numel(c), plot(c(k).dU,c(k).C,'Color',[0.6 0.6 0.6 0.5],'LineStyle','-','Marker','none','LineWidth',0.5,'HandleVisibility','off'); end
    for g = 1:3
        ok = ~isnan(G{:,['Gen' num2str(g) '_V']});
        scatter(G{ok,['Gen' num2str(g) '_V']},G{ok,['Gen' num2str(g) '_pF']},22,EvCol(g+1,:),'filled')
    end
    xlim(xl); xlabel('Swing since the previous extremum / V'); ylabel('dQ_p/dU_r / pF'); box on; grid on
    Leg = EvName(2:4);
    if Img.Ok
        yyaxis right
        plot(Img.Swing,Img.Tot,'k-','LineWidth',1.6)
        plot(Img.PkTot,interp1(Img.Swing,Img.Tot,Img.PkTot),'k^','MarkerFaceColor','k','MarkerSize',9)
        ylabel('ICCD 337 intensity, 3 kHz (normalised)'); ylim([0 1.6]); set(gca,'YColor','k')
        yyaxis left
        Leg = [Leg {'ICCD 337 intensity','ICCD intensity peaks'}];
    end
    legend(Leg,'Location','northwest')
    title('Peaks accepted inside the dQ_p/dU_r windows (shaded)')
    SaveFigure(fig,OutDir,['GenerationStaircase_' strrep(Sweeps(s).Name,' ','')])
end

%% Local functions
function [b,ci] = SlopeCI(x,y)
% Least-squares slope of y on x with its 95% confidence interval
ok = isfinite(x) & isfinite(y);  x = x(ok);  y = y(ok);  n = numel(x);
X = [ones(n,1) x];  c = X\y;  r = y - X*c;
se = sqrt(sum(r.^2)/(n-2)/sum((x - mean(x)).^2));
b = c(2);  ci = b + [-1 1]*tinv(0.975,n-2)*se;
end

function [Resr,Coff] = OffMeasurement(Folder,AmpGain,CurrentGain,f0)
% Mean equivalent series resistance and off capacitance over the off traces
% recorded within 5% of the nominal frequency
Files = dir(fullfile(Folder,'*.csv'));
R = []; C = [];
for k = 1:numel(Files)
    [t,V,I] = ReadScope(fullfile(Files(k).folder,Files(k).name),AmpGain,CurrentGain);
    f = FitFrequency(t,V,f0);
    if abs(f - f0) > 0.05*f0, continue, end
    Z = Phasor(t,V,f)/Phasor(t,I,f);
    R(end+1) = real(Z);                     %#ok<AGROW>
    C(end+1) = -1/(2*pi*f*imag(Z));         %#ok<AGROW>
end
if isempty(R), error('No off traces near %g Hz in %s.',f0,Folder); end
Resr = mean(R);  Coff = mean(C);
end

function [t,V,I] = ReadScope(File,AmpGain,CurrentGain)
% Tektronix MSO64 _ALL.csv: TIME, CH1 (voltage monitor), CH2 (current monitor)
M = readmatrix(File,'NumHeaderLines',11);
t = M(:,1);  V = M(:,2)*AmpGain;  I = M(:,3)*CurrentGain;
end

function f = FitFrequency(t,y,f0)
% Frequency that best fits a sine to y, starting from the nominal f0
B = @(f) [sin(2*pi*f*t) cos(2*pi*f*t) ones(size(t))];
f = fminsearch(@(f) norm(y - B(f)*(B(f)\y)),f0,optimset('TolX',1e-3,'Display','off'));
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
