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
OutDir   = fullfile(RepoRoot,'Data','FullImagingDataSet','Results','EoN_BigHex');
DelayThreshold = 4;          % us after the median side
ThresholdSweep = 2:2:16;     % us, to check the result does not hinge on the threshold
DelayBinEdges  = [-Inf 0 4 8 12 16 24 Inf];
NBoot          = 5000;       % bootstrap resamples for the median difference

D = load(fullfile(OutDir,'BigHexDerived.mat'),'SideTable');
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
