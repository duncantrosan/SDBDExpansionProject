%% Delayed vs prompt sides - pooled over the whole electrode
%
% Puts every hexagon side in one pot (no pairing within a hexagon) and
% compares the E/N of delayed sides with all the other sides, per half-wave,
% in the three windows from ICCDFilterRatioEoN.m Section 10:
%   Matched  - the same gates for every side, after 95% of sides have ignited
%   Onset    - the first gates after each side's own ignition
%   HalfWave - the whole half-wave
%
% Reads only the per-side table cached in BigHexDerived.mat, so the delay
% threshold can be changed and rerun in seconds.

clear
close all

%% Section 1 - Options
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
OutDir   = fullfile(RepoRoot,'Output','SideAnalysis');
if ~exist(OutDir,'dir'), mkdir(OutDir); end
DelayThreshold = 4;          % us after the median side
ThresholdSweep = 2:2:16;     % us, to check the result does not hinge on the threshold
DelayBinEdges  = [-Inf 0 4 8 12 16 24 Inf];
NBoot          = 5000;       % bootstrap resamples for the median difference

MinDepthCounts = 200;        % mean 337 counts per pixel per gate to show a depth bin
NBootDepth     = 1000;       % bootstrap resamples (over sides) for the depth profiles

D = load(fullfile(RepoRoot,'Output','Cache','BigHexDerived.mat'),'SideTable','SideProf','SideProfFNS','DepthMM', ...
         'MatchStart','HumpStart','HumpEnd','HalfWave','GateDelay','FilterFactor337', ...
         'FilterFactor391','EoNCoeff');
T = D.SideTable;
HWNames = {'Negative','Positive'};
Win     = {'Matched','Onset','HalfWave'};

%% Section 2 - Pooled comparison
Summary = table;
rng(1)
for hw = 1:2
    t = T(strcmp(T.HalfWave,HWNames{hw}),:);
    Ign = ~isnan(t.IgnitionTime_us);
    Delayed = Ign & t.Delay_us >= DelayThreshold;
    Prompt  = Ign & ~Delayed;
    fprintf('\n%s half-wave: %d sides ignite, %d delayed (>= %g us after the median side)\n', ...
            HWNames{hw},nnz(Ign),nnz(Delayed),DelayThreshold);
    for w = 1:3
        E = t.(['EoN' Win{w} '_Td']);
        p = E(Prompt & ~isnan(E));  d = E(Delayed & ~isnan(E));
        Boot = zeros(NBoot,1);
        for b = 1:NBoot
            Boot(b) = median(d(randi(numel(d),numel(d),1))) - median(p(randi(numel(p),numel(p),1)));
        end
        CI = prctile(Boot,[2.5 97.5]);
        Cliff = mean(sign(d(:) - p(:)'),'all');       % P(delayed > prompt) - P(delayed < prompt)
        pr = ranksum(d,p);
        fprintf('  %-9s prompt %5.0f Td | delayed %5.0f Td | diff %+5.0f Td [%+.0f, %+.0f] | Cliff''s delta %+.2f | ranksum p = %.2g\n', ...
                Win{w},median(p),median(d),median(d)-median(p),CI(1),CI(2),Cliff,pr);
        Summary = [Summary; table(HWNames(hw),Win(w),numel(p),numel(d),median(p),median(d), ...
                   median(d)-median(p),CI(1),CI(2),Cliff,pr,'VariableNames', ...
                   {'HalfWave','Window','NPrompt','NDelayed','EoNPrompt_Td','EoNDelayed_Td', ...
                    'Diff_Td','DiffLo_Td','DiffHi_Td','CliffsDelta','RanksumP'})]; %#ok<AGROW>
    end
    % Does the answer depend on where the delay threshold is put?
    fprintf('  Matched-window difference vs threshold: ');
    for thr = ThresholdSweep
        dl = Ign & t.Delay_us >= thr;
        fprintf('%g us: %+.0f Td (n=%d)  ',thr, ...
                median(t.EoNMatched_Td(dl),'omitnan') - median(t.EoNMatched_Td(Ign & ~dl),'omitnan'),nnz(dl));
    end
    fprintf('\n');
end
writetable(Summary,fullfile(OutDir,'SideDelayPooled.csv'));

%% Section 3 - Figure: one row per half-wave
fig = figure('Color','w','Position',[50 50 1500 820]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
Col = [0.3 0.5 0.9; 0.9 0.2 0.2];
for hw = 1:2
    t = T(strcmp(T.HalfWave,HWNames{hw}),:);
    Ign = ~isnan(t.IgnitionTime_us);
    Delayed = Ign & t.Delay_us >= DelayThreshold;

    % Box plots, all sides pooled, for each window
    nexttile
    X = []; Y = []; G = [];
    for w = 1:3
        E = t.(['EoN' Win{w} '_Td']);
        k = Ign & ~isnan(E);
        X = [X; w*ones(nnz(k),1)];  Y = [Y; E(k)];  G = [G; Delayed(k)]; %#ok<AGROW>
    end
    boxchart(categorical(X,1:3,Win),Y,'GroupByColor',categorical(G,[0 1],{'Prompt','Delayed'}));
    colororder(gca,Col); ylabel('E/N / Td'); legend('Location','best')
    title(sprintf('%s: all sides pooled (delayed >= %g \\mus)',HWNames{hw},DelayThreshold))

    % E/N against delay, binned
    nexttile
    Bin = discretize(t.Delay_us,DelayBinEdges);
    Mid = [-4 2 6 10 14 20 30];
    m = accumarray(Bin(Ign),t.EoNMatched_Td(Ign),[numel(Mid) 1],@(x) mean(x,'omitnan'),NaN);
    s = accumarray(Bin(Ign),t.EoNMatched_Td(Ign),[numel(Mid) 1],@(x) std(x,'omitnan')/sqrt(numel(x)),NaN);
    n = accumarray(Bin(Ign),1,[numel(Mid) 1]);
    errorbar(Mid,m,s,'ko-','LineWidth',1.3,'MarkerFaceColor','k'); hold on
    for b = find(n' > 0)
        text(Mid(b),m(b)+s(b)+3,sprintf('n=%d',n(b)),'HorizontalAlignment','center','FontSize',8)
    end
    xline(DelayThreshold,'--'); box on
    xlabel('Side ignition delay / \mus (after the median side)'); ylabel('Matched-window E/N / Td')
    title('E/N vs delay (mean \pm s.e.m.)')

    % Distributions
    nexttile
    [F1,x1] = ecdf(t.EoNMatched_Td(Ign & ~Delayed));
    [F2,x2] = ecdf(t.EoNMatched_Td(Delayed));
    stairs(x1,F1,'Color',Col(1,:),'LineWidth',1.5); hold on
    stairs(x2,F2,'Color',Col(2,:),'LineWidth',1.5); box on
    xlabel('Matched-window E/N / Td'); ylabel('Cumulative fraction of sides')
    legend('Prompt','Delayed','Location','southeast'); title('Distribution, all sides')
end
exportgraphics(fig,fullfile(OutDir,'SideDelayPooled.png'),'Resolution',150)

%% Section 4 - E/N against distance from the edge
% The apparent E/N of a side is weighted by where it emits. Delayed sides emit
% relatively less right at the edge, where E/N is highest, so part of their
% lower apparent E/N could be weighting rather than a lower local field. Here
% the ratio is formed at each depth bin (summed over the matched gates and over
% the sides in each group), and the gap is split in two:
%   R = sum_d w_d r_d   (w_d = share of the 337 light at depth d, r_d = ratio there)
%   weighting effect   = E/N( sum_d wDelayed_d * rPrompt_d ) - E/N( RPrompt )
%   local-field effect = E/N( RDelayed ) - E/N( sum_d wDelayed_d * rPrompt_d )
DepthRatio = @(A,F) (sum(F,1)/D.FilterFactor391)./(sum(A,1)/D.FilterFactor337);
EoN = @(r) polyval(D.EoNCoeff,r);
fig2 = figure('Color','w','Position',[60 60 1300 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
Depth = D.DepthMM(:)';
DepthSummary = table;
rng(2)
for hw = 1:2
    h  = find((D.HalfWave(D.HumpStart) > 0) == (hw == 2),1);
    gm = D.MatchStart(h):D.HumpEnd(h);
    t  = T(strcmp(T.HalfWave,HWNames{hw}),:);
    Ign = ~isnan(t.IgnitionTime_us);
    Grp = {Ign & t.Delay_us < DelayThreshold, Ign & t.Delay_us >= DelayThreshold};
    A = sum(D.SideProf(:,:,gm,2),3);          % sides x depth, 337 summed over matched gates
    F = sum(D.SideProfFNS(:,:,gm),3);         % same for the corrected N2+ signal

    r = zeros(2,numel(Depth));  w = r;  rLo = r;  rHi = r;
    Show = true(1,numel(Depth));
    for gi = 1:2
        Ag = A(Grp{gi},:);  Fg = F(Grp{gi},:);
        r(gi,:) = DepthRatio(Ag,Fg);
        w(gi,:) = sum(Ag,1);
        Show = Show & sum(Ag,1)/nnz(Grp{gi})/numel(gm) >= MinDepthCounts;
        Boot = zeros(NBootDepth,numel(Depth));
        idx = find(Grp{gi});
        for b = 1:NBootDepth
            s = idx(randi(numel(idx),numel(idx),1));
            Boot(b,:) = EoN(DepthRatio(A(s,:),F(s,:)));
        end
        rLo(gi,:) = prctile(Boot,2.5);  rHi(gi,:) = prctile(Boot,97.5);
    end
    w = w(:,Show)./sum(w(:,Show),2);          % emission share over the shown bins
    rs = r(:,Show);
    RP = sum(w(1,:).*rs(1,:));  RD = sum(w(2,:).*rs(2,:));  RX = sum(w(2,:).*rs(1,:));
    Weighting = EoN(RX) - EoN(RP);
    LocalField = EoN(RD) - EoN(RX);
    fprintf('\n%s half-wave, depth-resolved (matched gates %d-%d us, %.2f to %.2f mm from the edge):\n', ...
            HWNames{hw},D.GateDelay(gm(1)),D.GateDelay(gm(end)),min(Depth(Show)),max(Depth(Show)));
    fprintf('  Apparent E/N: prompt %.0f Td, delayed %.0f Td, gap %+.0f Td\n',EoN(RP),EoN(RD),EoN(RD)-EoN(RP));
    fprintf('    from emission weighting : %+.0f Td\n',Weighting);
    fprintf('    from local field        : %+.0f Td\n',LocalField);
    for dd = [0 0.25 0.5 0.75 1.0 1.25]
        [~,j] = min(abs(Depth - dd));
        if Show(j)
            fprintf('  at %.2f mm: prompt %.0f Td [%.0f-%.0f], delayed %.0f Td [%.0f-%.0f]\n',Depth(j), ...
                    EoN(r(1,j)),rLo(1,j),rHi(1,j),EoN(r(2,j)),rLo(2,j),rHi(2,j));
        end
    end
    DepthSummary = [DepthSummary; table(HWNames(hw),EoN(RP),EoN(RD),Weighting,LocalField, ...
        'VariableNames',{'HalfWave','EoNPrompt_Td','EoNDelayed_Td','FromWeighting_Td','FromLocalField_Td'})]; %#ok<AGROW>

    nexttile
    Col2 = [0.3 0.5 0.9; 0.9 0.2 0.2];
    for gi = 1:2
        x = Depth(Show);
        fill([x fliplr(x)],[rLo(gi,Show) fliplr(rHi(gi,Show))],Col2(gi,:),'FaceAlpha',0.15, ...
             'EdgeColor','none','HandleVisibility','off'); hold on
        plot(x,EoN(r(gi,Show)),'o-','Color',Col2(gi,:),'MarkerFaceColor',Col2(gi,:),'MarkerSize',4,'LineWidth',1.3)
    end
    xline(0,'k:','Foil edge','HandleVisibility','off'); box on
    xlabel('Distance from the electrode edge / mm'); ylabel('E/N / Td')
    legend(sprintf('Prompt (n=%d)',nnz(Grp{1})),sprintf('Delayed (n=%d)',nnz(Grp{2})),'Location','best')
    title(sprintf('%s: E/N vs depth, matched gates (95%% CI over sides)',HWNames{hw}))
    nexttile
    plot(Depth(Show),w(1,:),'o-','Color',Col2(1,:),'MarkerSize',4,'LineWidth',1.3); hold on
    plot(Depth(Show),w(2,:),'o-','Color',Col2(2,:),'MarkerSize',4,'LineWidth',1.3)
    xline(0,'k:','HandleVisibility','off'); box on
    xlabel('Distance from the electrode edge / mm'); ylabel('Share of the 337 light')
    legend('Prompt','Delayed','Location','best')
    title(sprintf('Where the light comes from: gap %+.0f Td = %+.0f weighting %+.0f local field', ...
          EoN(RD)-EoN(RP),Weighting,LocalField))
end
exportgraphics(fig2,fullfile(OutDir,'SideDelayDepth.png'),'Resolution',150)
writetable(DepthSummary,fullfile(OutDir,'SideDelayDepth.csv'));
