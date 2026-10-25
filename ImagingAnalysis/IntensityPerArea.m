%% Intensity per area - brighter emission or more plasma?
%
% Splits the total 337 intensity of each half-wave into the plasma area and
% the intensity per area,
%   I = A * J,   ln(I/I0) = ln(A/A0) + ln(J/J0),
% with I the 337 counts inside the plasma mask, A the mask area and J = I/A the
% mean 337 counts per cm^2 of lit plasma. The log split says, from the
% reference gate on, how much of the growth of I is expansion (A) and how much
% is stronger emission (J).
%
% Reads only BigHexDerived.mat (ICCDFilterRatioEoN.m, Section 10). The mask
% uses the PerimeterCodev2.m thresholds, which differ by half-wave (3000 counts
% negative, 1580 positive), so compare J within a half-wave, not across.
% A thresholded mask also grows when existing plasma brightens past the
% threshold, so part of what lands in A at low intensity is really J.
%
% Positions are the swing since the previous Ur extremum, as in
% BreakdownGenerationTest.m, with the image intensity peaks marked.

clear
close all

%% Section 1 - Options
RepoRoot = fileparts(fileparts(mfilename('fullpath')));
OutDir   = fullfile(RepoRoot,'Output','IntensityPerArea');
if ~exist(OutDir,'dir'), mkdir(OutDir); end
D = load(fullfile(RepoRoot,'Output','Cache','BigHexDerived.mat'), ...
         'UrGate','HalfWave','HumpStart','HumpEnd','GateDelay','PlasmaArea','TotalInt','Tot337','IgnitionCounts');

MinAreaFrac   = 0.02;   % J only where A is above this share of the half-wave maximum
PeakProminence = 0.03;  % of the normalised intensity, as in BreakdownGenerationTest.m

%% Section 2 - Per half-wave: I, A, J and the log split
HWNames = {'Negative','Positive'};
Out = table;
for h = 1:numel(D.HumpStart)
    g  = D.HumpStart(h):D.HumpEnd(h);
    hw = 1 + (D.HalfWave(g(1)) > 0);
    if hw == 2, Swing = D.UrGate(g) - min(D.UrGate); else, Swing = max(D.UrGate) - D.UrGate(g); end
    I = D.TotalInt(g);
    A = D.PlasmaArea(g)/100;                 % cm^2
    J = I./A;
    J(A < MinAreaFrac*max(A)) = NaN;
    r = find(~isnan(J),1);                   % reference gate: first with a usable J
    LnI = log(I/I(r));  LnA = log(A/A(r));  LnJ = log(J/J(r));
    [~,lp] = findpeaks(I/max(I),'MinPeakProminence',PeakProminence);

    H(hw) = struct('Swing',Swing,'I',I,'A',A,'J',J,'LnI',LnI,'LnA',LnA,'LnJ',LnJ,'Ref',r,'Pk',lp); %#ok<SAGROW>
    Out = [Out; table(repmat(string(HWNames{hw}),numel(g),1),D.GateDelay(g)',D.UrGate(g)',Swing', ...
           I',A',J',LnI',LnA',LnJ','VariableNames',{'HalfWave','GateDelay_us','Ur_V','Swing_V', ...
           'TotalInt337','Area_cm2','IntPerArea_cm2','LnI','LnA','LnJ'})]; %#ok<AGROW>

    fprintf('\n%s half-wave (mask threshold %d counts), reference gate at swing %.0f V:\n', ...
            HWNames{hw},D.IgnitionCounts(hw),Swing(r));
    [~,iA] = max(A);  [~,iJ] = max(J);
    fprintf('  area max %.2f cm^2 at %.0f V, intensity/area max %.3g counts/cm^2 at %.0f V\n', ...
            A(iA),Swing(iA),J(iJ),Swing(iJ));
    for k = lp
        fprintf('  intensity peak at %4.0f V: I x%6.1f = A x%5.1f * J x%4.2f  (expansion %3.0f%%, emission %3.0f%% of ln growth)\n', ...
                Swing(k),exp(LnI(k)),exp(LnA(k)),exp(LnJ(k)),100*LnA(k)/LnI(k),100*LnJ(k)/LnI(k));
    end
end
writetable(Out,fullfile(OutDir,'IntensityPerArea.csv'));

%% Section 3 - Figure: rows I & A, J, log split; columns negative, positive
fig = figure('Color','w','Position',[50 40 1300 1050]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
Col = [0.15 0.3 0.8; 0.85 0.3 0.1];
for row = 1:3
    for hw = 1:2
        S = H(hw);  x = S.Swing;
        nexttile; hold on
        switch row
            case 1
                yyaxis left
                plot(x,S.I/1e9,'o-','Color',Col(hw,:),'MarkerSize',3,'LineWidth',1.4)
                ylabel('Total 337 intensity in mask / 10^9 counts'); set(gca,'YColor',Col(hw,:))
                yyaxis right
                plot(x,S.A,'s-','Color',[0.2 0.6 0.3],'MarkerSize',3,'LineWidth',1.4)
                ylabel('Plasma area / cm^2'); set(gca,'YColor',[0.2 0.6 0.3])
                title(sprintf('%s half-wave',HWNames{hw}))
                legend({'Intensity','Area'},'Location','northwest')
            case 2
                plot(x,S.J/1e7,'o-','Color','k','MarkerSize',3,'LineWidth',1.4)
                ylabel('Intensity / area / 10^7 counts cm^{-2}')
                ylim([0 max(S.J/1e7,[],'omitnan')*1.1])
            case 3
                ok = ~isnan(S.LnJ);
                plot(x(ok),S.LnI(ok),'-','Color',Col(hw,:),'LineWidth',1.8)
                plot(x(ok),S.LnA(ok),'-','Color',[0.2 0.6 0.3],'LineWidth',1.4)
                plot(x(ok),S.LnJ(ok),'-','Color','k','LineWidth',1.4)
                yline(0,'Color',[0.6 0.6 0.6],'HandleVisibility','off')
                ylabel(sprintf('ln(X / X at %.0f V)',x(S.Ref)))
                legend({'ln I (total)','ln A (expansion)','ln J (intensity / area)'},'Location','south')
        end
        for k = S.Pk, xline(x(k),':','Color',[0.4 0.4 0.4],'HandleVisibility','off'); end
        xlabel('Swing since the previous extremum / V'); box on; grid on
        xlim([min(x) max(x)])
        ax = gca;  ax.Toolbar.Visible = 'off';
        text(0.98,0.92,char('A' + 2*(row-1) + hw - 1),'Units','normalized','FontWeight','bold', ...
             'FontSize',14,'HorizontalAlignment','right')
    end
end
fig.Name = 'IntensityPerArea';
exportgraphics(fig,fullfile(OutDir,'IntensityPerArea.png'),'Resolution',200)
savefig(fig,fullfile(OutDir,'IntensityPerArea.fig'))
fprintf('\nFigure and Origin table written to %s\n',OutDir);
