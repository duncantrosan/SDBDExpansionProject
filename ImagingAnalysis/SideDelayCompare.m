%% Delayed vs prompt sides - pooled over the whole electrode
%
% Puts every hexagon side in one pot (no pairing within a hexagon) and
% compares the N2(C)/N2+(B) emission ratio of delayed sides with all the other
% sides, per half-wave, in the three windows from ICCDFilterRatioEoN.m
% Section 10:
%   Matched  - the same gates for every side, after 95% of sides have ignited
%   Onset    - the first gates after each side's own ignition
%   HalfWave - the whole half-wave
% The ratio is echelle-anchored (see ICCDFilterRatioEoN.m). A higher
% N2(C)/N2+(B) means a lower field.
%
% Reads only the per-side data cached in BigHexDerived.mat, so the delay
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
NBoot          = 5000;       % bootstrap resamples for the ratio of medians
ReportEoN      = 0;          % 1 = also print E/N (Td) from the calibration table

MinDepthCounts = 200;        % mean 337 counts per pixel per gate to show a depth bin
NBootDepth     = 1000;       % bootstrap resamples (over sides) for the depth profiles

D = load(fullfile(RepoRoot,'Output','Cache','BigHexDerived.mat'),'SideTable','SideProf','SideProfFNS','DepthMM', ...
         'MatchStart','HumpStart','HumpEnd','HalfWave','GateDelay','FilterFactor337', ...
         'FilterFactor391','EoNCoeff');
T = D.SideTable;
HWNames = {'Negative','Positive'};
Win     = {'Matched','Onset','HalfWave'};
YLab    = 'N_2(C) / N_2^+(B)';
EoN     = @(cb) polyval(D.EoNCoeff,1./cb);

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
        V = t.(['CB' Win{w}]);
        p = V(Prompt & ~isnan(V));  d = V(Delayed & ~isnan(V));
        Boot = zeros(NBoot,1);
        for b = 1:NBoot
            Boot(b) = median(d(randi(numel(d),numel(d),1)))/median(p(randi(numel(p),numel(p),1))) - 1;
        end
        CI = 100*prctile(Boot,[2.5 97.5]);
        Diff = 100*(median(d)/median(p) - 1);
        Cliff = mean(sign(d(:) - p(:)'),'all');       % P(delayed > prompt) - P(delayed < prompt)
        pr = ranksum(d,p);
        fprintf('  %-9s N2(C)/N2+(B): prompt %.2f | delayed %.2f | %+.1f%% [%+.1f, %+.1f] | Cliff''s delta %+.2f | ranksum p = %.2g', ...
                Win{w},median(p),median(d),Diff,CI(1),CI(2),Cliff,pr);
        if ReportEoN == 1
            fprintf(' | E/N %.0f vs %.0f Td',EoN(median(p)),EoN(median(d)));
        end
        fprintf('\n');
        Summary = [Summary; table(HWNames(hw),Win(w),numel(p),numel(d),median(p),median(d), ...
                   Diff,CI(1),CI(2),Cliff,pr,'VariableNames', ...
                   {'HalfWave','Window','NPrompt','NDelayed','CBPrompt','CBDelayed', ...
                    'Diff_pct','DiffLo_pct','DiffHi_pct','CliffsDelta','RanksumP'})]; %#ok<AGROW>
    end
    % Does the answer depend on where the delay threshold is put?
    fprintf('  Matched-window difference vs threshold: ');
    for thr = ThresholdSweep
        dl = Ign & t.Delay_us >= thr;
        fprintf('%g us: %+.1f%% (n=%d)  ',thr, ...
                100*(median(t.CBMatched(dl),'omitnan')/median(t.CBMatched(Ign & ~dl),'omitnan') - 1),nnz(dl));
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
        V = t.(['CB' Win{w}]);
        k = Ign & ~isnan(V);
        X = [X; w*ones(nnz(k),1)];  Y = [Y; V(k)];  G = [G; Delayed(k)]; %#ok<AGROW>
    end
    boxchart(categorical(X,1:3,Win),Y,'GroupByColor',categorical(G,[0 1],{'Prompt','Delayed'}));
    colororder(gca,Col); ylabel(YLab); legend('Location','best')
    title(sprintf('%s: all sides pooled (delayed >= %g \\mus)',HWNames{hw},DelayThreshold))

    % Ratio against delay, binned
    nexttile
    Bin = discretize(t.Delay_us,DelayBinEdges);
    Mid = [-4 2 6 10 14 20 30];
    m = accumarray(Bin(Ign),t.CBMatched(Ign),[numel(Mid) 1],@(x) mean(x,'omitnan'),NaN);
    s = accumarray(Bin(Ign),t.CBMatched(Ign),[numel(Mid) 1],@(x) std(x,'omitnan')/sqrt(numel(x)),NaN);
    n = accumarray(Bin(Ign),1,[numel(Mid) 1]);
    errorbar(Mid,m,s,'ko-','LineWidth',1.3,'MarkerFaceColor','k'); hold on
    for b = find(n' > 0)
        text(Mid(b),m(b)+1.15*s(b),sprintf('n=%d',n(b)),'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom','FontSize',8)
    end
    xline(DelayThreshold,'--'); box on
    xlabel('Side ignition delay / \mus (after the median side)'); ylabel(['Matched-window ' YLab])
    title('Ratio vs delay (mean \pm s.e.m.)')

    % Distributions
    nexttile
    [F1,x1] = ecdf(t.CBMatched(Ign & ~Delayed));
    [F2,x2] = ecdf(t.CBMatched(Delayed));
    stairs(x1,F1,'Color',Col(1,:),'LineWidth',1.5); hold on
    stairs(x2,F2,'Color',Col(2,:),'LineWidth',1.5); box on
    xlabel(['Matched-window ' YLab]); ylabel('Cumulative fraction of sides')
    legend('Prompt','Delayed','Location','southeast'); title('Distribution, all sides')
end
exportgraphics(fig,fullfile(OutDir,'SideDelayPooled.png'),'Resolution',150)

%% Section 4 - Ratio against distance from the edge
% The apparent ratio of a side is weighted by where it emits. Delayed sides
% emit relatively less right at the edge, where the field is highest, so part
% of their higher apparent N2(C)/N2+(B) could be weighting rather than a lower
% local field. Here the ratio is formed at each depth bin (summed over the
% matched gates and over the sides in each group), and the gap is split in two.
% Working in N2+(B)/N2(C), which adds up linearly with the 337 light:
%   R = sum_d w_d r_d   (w_d = share of the 337 light at depth d, r_d = N2+/N2 there)
%   gap in N2(C)/N2+(B)  = RPrompt/RDelayed - 1
%     = (1 + weighting) * (1 + local field) - 1, with
%   weighting   = RPrompt / sum_d wDelayed_d*rPrompt_d - 1
%   local field = sum_d wDelayed_d*rPrompt_d / RDelayed - 1
DepthNB = @(A,F) (sum(F,1)/D.FilterFactor391)./(sum(A,1)/D.FilterFactor337);   % N2+(B)/N2(C)
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
    F = sum(D.SideProfFNS(:,:,gm),3);         % same for the corrected, anchored N2+ signal

    r = zeros(2,numel(Depth));  w = r;  cLo = r;  cHi = r;
    Show = true(1,numel(Depth));
    for gi = 1:2
        Ag = A(Grp{gi},:);  Fg = F(Grp{gi},:);
        r(gi,:) = DepthNB(Ag,Fg);
        w(gi,:) = sum(Ag,1);
        Show = Show & sum(Ag,1)/nnz(Grp{gi})/numel(gm) >= MinDepthCounts;
        Boot = zeros(NBootDepth,numel(Depth));
        idx = find(Grp{gi});
        for b = 1:NBootDepth
            s = idx(randi(numel(idx),numel(idx),1));
            Boot(b,:) = 1./DepthNB(A(s,:),F(s,:));
        end
        cLo(gi,:) = prctile(Boot,2.5);  cHi(gi,:) = prctile(Boot,97.5);
    end
    CBd = 1./r;
    w = w(:,Show)./sum(w(:,Show),2);          % emission share over the shown bins
    rs = r(:,Show);
    RP = sum(w(1,:).*rs(1,:));  RD = sum(w(2,:).*rs(2,:));  RX = sum(w(2,:).*rs(1,:));
    Gap        = 100*(RP/RD - 1);
    Weighting  = 100*(RP/RX - 1);
    LocalField = 100*(RX/RD - 1);
    fprintf('\n%s half-wave, depth-resolved (matched gates %d-%d us, %.2f to %.2f mm from the edge):\n', ...
            HWNames{hw},D.GateDelay(gm(1)),D.GateDelay(gm(end)),min(Depth(Show)),max(Depth(Show)));
    fprintf('  Apparent N2(C)/N2+(B): prompt %.2f, delayed %.2f, gap %+.1f%%\n',1/RP,1/RD,Gap);
    fprintf('    from emission weighting : %+.1f%%\n',Weighting);
    fprintf('    from local field        : %+.1f%%\n',LocalField);
    for dd = [0 0.25 0.5 0.75 1.0 1.25]
        [~,j] = min(abs(Depth - dd));
        if Show(j)
            fprintf('  at %.2f mm: prompt %.2f [%.2f-%.2f], delayed %.2f [%.2f-%.2f]\n',Depth(j), ...
                    CBd(1,j),cLo(1,j),cHi(1,j),CBd(2,j),cLo(2,j),cHi(2,j));
        end
    end
    DepthSummary = [DepthSummary; table(HWNames(hw),1/RP,1/RD,Gap,Weighting,LocalField, ...
        'VariableNames',{'HalfWave','CBPrompt','CBDelayed','Gap_pct','FromWeighting_pct','FromLocalField_pct'})]; %#ok<AGROW>

    nexttile
    Col2 = [0.3 0.5 0.9; 0.9 0.2 0.2];
    for gi = 1:2
        x = Depth(Show);
        fill([x fliplr(x)],[cLo(gi,Show) fliplr(cHi(gi,Show))],Col2(gi,:),'FaceAlpha',0.15, ...
             'EdgeColor','none','HandleVisibility','off'); hold on
        plot(x,CBd(gi,Show),'o-','Color',Col2(gi,:),'MarkerFaceColor',Col2(gi,:),'MarkerSize',4,'LineWidth',1.3)
    end
    xline(0,'k:','Foil edge','HandleVisibility','off'); box on
    xlabel('Distance from the electrode edge / mm'); ylabel(YLab)
    legend(sprintf('Prompt (n=%d)',nnz(Grp{1})),sprintf('Delayed (n=%d)',nnz(Grp{2})),'Location','best')
    title(sprintf('%s: ratio vs depth, matched gates (95%% CI over sides)',HWNames{hw}))
    nexttile
    plot(Depth(Show),w(1,:),'o-','Color',Col2(1,:),'MarkerSize',4,'LineWidth',1.3); hold on
    plot(Depth(Show),w(2,:),'o-','Color',Col2(2,:),'MarkerSize',4,'LineWidth',1.3)
    xline(0,'k:','HandleVisibility','off'); box on
    xlabel('Distance from the electrode edge / mm'); ylabel('Share of the 337 light')
    legend('Prompt','Delayed','Location','best')
    title(sprintf('Gap %+.1f%% = %+.1f%% weighting, %+.1f%% local field',Gap,Weighting,LocalField))
end
exportgraphics(fig2,fullfile(OutDir,'SideDelayDepth.png'),'Resolution',150)
writetable(DepthSummary,fullfile(OutDir,'SideDelayDepth.csv'));
