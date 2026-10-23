%% Surface Dielectric Barrier Discharge Analysis
%
% Interprets the electrical characteristics of a surface dielectric barrier
% discharge (SDBD) and returns plasma parameters: electron density, reduced
% electric field, discharge current, dielectric capacitance and gap voltage.
% The plasma is modeled with the equivalent circuit of Figure 1, a modified
% form of the classic volume-DBD circuit (Manley, 1943).
%
% DEPENDENCIES
%   OscopeReadin        - reads Tektronix oscilloscope CSV exports
%                         (developed for the MDO31604; should work for most
%                         Tektronix scopes saving to CSV)
%   Curve Fitting Toolbox     - fit, fittype, fitoptions, prepareCurveData,
%                               coeffvalues, smooth, spaps
%   Signal Processing Toolbox - sgolay, downsample
%   MATLAB R2018b or newer    - xline
%
% EXPECTED INPUTS
%   PathOff - directory of Tektronix CSVs taken with the plasma OFF (low
%             voltage). A directory is required, not a single file: the code
%             averages several traces to get the equivalent series resistance
%             and off capacitance plus a standard deviation for error
%             propagation. Three files minimum is recommended.
%   PathOn  - directory of CSVs with the plasma ON. Traces below ignition are
%             fine here and are often useful for locating the ignition point.
%   Optional: a reduced electric field from an independent measurement (e.g.
%             optical emission spectroscopy) or a secondary electron emission
%             coefficient. One of these is needed to build a Paschen curve and
%             locate the burning voltage on it. This step can be skipped, but
%             reduced electric field and electron density are then unavailable.
%
% OUTPUTS (left in the workspace)
%   R_off / ResistanceAverage   - equivalent series resistance, Ohms
%   C_off / CapacitanceAverage  - off capacitance, F
%   Discharge                   - discharge current, A
%   DielectricCapacitancePositive / Negative - F
%   PlasmaResistancePositive / Negative      - Ohms
%   UgTrunc                     - gap voltage vs. time
%   ChargeDerivitive            - dQ/dV vs. time, F
%
% LIMITATIONS
%   Assumes a pure sine-wave excitation. Multiple driving frequencies and
%   complicated geometries are not yet supported.
%
% ---------------------------------------------------------------------------
% Converted from SurfaceDBDSimpleLinearModelV3.mlx to a plain .m script.
% Lines marked % FIX: were changed from the original. Lines marked % REVIEW:
% are physics/intent questions left for you to decide - behavior unchanged.
% ---------------------------------------------------------------------------

%% Section 1 - Configuration and plasma-off read in
clear
close all

% --- Logic flags ---------------------------------------------------------
% Set for ease of use. Most exist to skip parts of the code or to hard code
% values instead of prompting.
OptionalInputs   = 1;   % 1 = dialogue box for measured gamma or reduced E field
VerifyParameters = 0;   % 1 = dialogue boxes to verify/change input parameters
SingleFrequency  = 1;   % 1 = single drive frequency (leave at 1 until multi-frequency is supported)
DeviceSelection  = 1;   % 0 = hard coded values, 1 = device selection dialogue, 2 = skip dependent results
BolsigResults    = '';  % Bolsig+ output used to determine electron density

% --- Numerical parameters ------------------------------------------------
SmoothWindow      = 41;        % FIX: was 40. sgolay needs an odd span; smooth()
                               % silently drops even spans by 1, so make it explicit.
ThreshInput       = 0.00035;   % A. FIX: comment said mA, but Discharge is in
                               % amps, so this is 0.35 mA.
Perimeter         = 0.2266;    % m,     perimeter of the powered electrode (Q/P equation)
CapacitancePerArea= 1.245e-7;  % F/m^2, needed for E field via the alpha equation
Thickness         = 90e-6;     % m,     dielectric thickness, gives E from V0

% --- Probe / amplifier calibration --------------------------------------
% FIX: pulled out of the read-in blocks so both On and Off use one definition.
AmpGain       = 500;      % V/V,  lab amplifier
CurrentGain   = 6/1000;   % A/V,  current monitor, 1 V = 6 mA

% --- Optional geometry for plasma extent ---------------------------------
ComputePlasmaExtent = 0;   % FIX: was commented out in the original but the
                           % trip-averaging loop at the end still referenced
                           % PlasmaExtentPositive, which threw an undefined
                           % variable error. Set to 1 and check Aside/Bside
                           % against your device before trusting the numbers.
Aside = 0.0896;            % m
Bside = 0.0237;            % m
n0    = 2.686e25;          % m^-3, Loschmidt number for Townsend conversion

% --- Paths ---------------------------------------------------------------
% REVIEW: the original mixed a Windows data path with macOS defaults for
% uigetdir and SaveFileName. Point all of these at one machine.
PathOff = 'E:\MATLAB\Data\SDBDModelingProject\VoltageandCurrent\2025\1600 Hz\bighex\bighex\vsweep\PlasmaOff';
PathOn  = 'E:\MATLAB\Data\SDBDModelingProject\VoltageandCurrent\2025\1600 Hz\bighex\bighex\vsweep\PlasmaOn';
DefaultBrowseDir = pwd;   % FIX: was a hard coded /Users/... path

RunName      = 'Hex Old';
SaveFileName = fullfile(pwd,'HexLattice.csv');   % FIX: added a semicolon and an
                                                 % extension. writecell needs one.

% --- Read the plasma-off traces ------------------------------------------
if isempty(PathOff)
    PathOff = uigetdir(DefaultBrowseDir,'Choose Plasma Off Directory');
end
[Time,Channel1,Channel2] = OscopeReadin(PathOff);
[Nfiles,Npoints] = size(Time);

if VerifyParameters == 1
    % FIX: the original left this branch empty, so setting VerifyParameters = 1
    % left Voltage and Current undefined and the next loop failed. Until the UI
    % exists, fall through to the calibrated values and warn.
    warning('VerifyParameters UI is not implemented; using the hard coded gains.');
end
Voltage = Channel1*AmpGain;
Current = Channel2*CurrentGain;

clear Channel1 Channel2

%% Section 2 - Off measurements: series resistance and off capacitance
%
% The off measurement is a simple series RC circuit (Figure 2). Fitting the
% voltage and current to sine waves and dividing their phasors gives a complex
% impedance; the real part is the equivalent series resistance and the
% imaginary part gives the off capacitance.
%
% Because both traces are reconstructed as ideal sinusoids, the impedance
% should be constant. A large standard deviation means the sine fit failed.
%
% NOTE on sign convention: vi and ii below are built as the complex conjugate
% of the usual phasor (exp(+1i*(pi/2 - theta)) rather than exp(1i*(theta - pi/2))).
% That flips the sign of imag(Z), which is why C_off = 1/(w*imag(Z)) rather
% than the textbook -1/(w*imag(Z)). The two sign errors cancel and C_off comes
% out positive, but do not "fix" one without the other.

% --- Preallocate ---------------------------------------------------------
F_meas1        = zeros(1,Nfiles);
F_meas2        = zeros(1,Nfiles);
VoltageMeasuredOff = zeros(1,Nfiles);   % FIX: renamed. The original reused
                                        % VoltageMeasured for both the off and
                                        % on runs, so if the on directory held
                                        % fewer files, stale off amplitudes
                                        % leaked into every plot in Section 7.
I_meas         = zeros(1,Nfiles);
phi1           = zeros(1,Nfiles);
phi2           = zeros(1,Nfiles);
Z_off          = zeros(1,Nfiles);
R_off          = zeros(1,Nfiles);
C_off          = zeros(1,Nfiles);
PhaseShiftTime = zeros(1,Nfiles);       % FIX: was growing inside the loop
indexOff1      = zeros(1,Nfiles);       % FIX: renamed, was shadowing the on-cycle index1
indexOff2      = zeros(1,Nfiles);
TruncIndexOff  = zeros(1,Nfiles);

ft   = fittype('sin1');
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';
opts.Lower   = [-Inf 0 -Inf];

for n = 1:Nfiles
    % First pass fit, used only to get an approximate frequency
    [xData,yData] = prepareCurveData(Time(n,:),Voltage(n,:));
    fitresult = fit(xData,yData,ft,opts);
    Vcoeff = coeffvalues(fitresult);
    F_meas1(n) = Vcoeff(2)/2/pi;

    % Index points for one full period (kept for duty-cycle work)
    [~,indexOff1(n)] = min(abs(Time(n,:) - 0/(4*F_meas1(n))));
    [~,indexOff2(n)] = min(abs(Time(n,:) - 4/(4*F_meas1(n))));
    TruncIndexOff(n) = abs(indexOff1(n) - indexOff2(n));

    % The off analysis uses the full record, not the truncated one
    TTruncOff = Time(n,:);
    VTruncOff = Voltage(n,:);
    ITruncOff = Current(n,:);

    % Refit the voltage
    [xData,yData] = prepareCurveData(TTruncOff,VTruncOff);
    fitresult = fit(xData,yData,ft,opts);
    Vcoeff = coeffvalues(fitresult);
    a1 = Vcoeff(1);  b1 = Vcoeff(2);  c1 = Vcoeff(3);

    % FIX: sin1 can return a negative amplitude with a phase 180 degrees off.
    % Normalise so the amplitude is positive, otherwise Z picks up a sign flip.
    if a1 < 0
        a1 = -a1;
        c1 = c1 + pi;
    end
    F_meas1(n)            = b1/2/pi;   % measured frequency, Hz
    VoltageMeasuredOff(n) = a1;        % voltage amplitude, V
    phi1(n)               = c1;        % phase

    % Fit the current
    [xData,yData] = prepareCurveData(TTruncOff,ITruncOff);
    fitresult = fit(xData,yData,ft,opts);
    Icoeff = coeffvalues(fitresult);
    a2 = Icoeff(1);  b2 = Icoeff(2);  c2 = Icoeff(3);
    if a2 < 0            % FIX: same amplitude-sign normalisation
        a2 = -a2;
        c2 = c2 + pi;
    end
    F_meas2(n) = b2/2/pi;
    I_meas(n)  = a2;      % current amplitude, A
    phi2(n)    = c2;
    PhaseShiftTime(n) = (1/F_meas1(n)/360)*rad2deg(abs(phi2(n) - phi1(n) - pi/2));

    % Rebuild the waveforms as phasors and take the impedance
    w  = F_meas1(n)*2*pi;
    vi = VoltageMeasuredOff(n)*exp(1i*(pi/2 - (w*Time(n,:) + phi1(n))));
    ii = I_meas(n)*exp(1i*(pi/2 - (w*Time(n,:) + phi2(n))));
    Z_off(n) = mean(vi./ii);

    err = std(vi./ii);
    % FIX: the original compared an absolute spread against 1e-7 regardless of
    % impedance magnitude. Made it relative so the check scales with Z.
    if err/abs(Z_off(n)) >= 1e-7
        fprintf('WARNING: std of Z is too high on file %d (relative spread %.2e)\n', ...
                n, err/abs(Z_off(n)));
    end

    R_off(n) = real(Z_off(n));
    C_off(n) = 1/(imag(Z_off(n))*w);   % see sign-convention note above

    if n ~= Nfiles
        clear VTruncOff ITruncOff TTruncOff
    end
end

R_off_err = std(R_off);
C_off_err = std(C_off);

figure
bar(R_off)
xlabel('Trial Number')
ylabel('Resistance / \Omega')
title('Series Resistance')

figure
bar(C_off*1e12)
xlabel('Trial Number')
ylabel('Capacitance / pF')
title('Off Capacitance')

fprintf('Average Equivalent Series Resistance = %.2f Ohms +- %.2f Ohms\n', ...
        mean(R_off), R_off_err);
fprintf('Average Off Capacitance = %.2f pF +- %.2f pF\n', ...
        mean(C_off)*1e12, C_off_err*1e12);

% Example: measured vs. simulated traces for the last file
figure
plot(TTruncOff,VTruncOff); hold on; plot(TTruncOff,real(vi))   % FIX: plot() on a
ylabel('Voltage / V')                                          % complex vector silently
yyaxis right                                                   % draws imag vs real.
plot(TTruncOff,ITruncOff); hold on; plot(TTruncOff,real(ii))   % Take the real part.
ylabel('Current / A')                                          % FIX: label was 'Current / I'
xlabel('Time / s')
legend('Measured Voltage','Simulated Voltage','Measured Current','Simulated Current')

ResistanceAverage  = mean(R_off);
CapacitanceAverage = mean(C_off);
ZZ = Z_off;

%% Section 3 - Plasma on read in
if isempty(PathOn)
    PathOn = uigetdir(DefaultBrowseDir,'Choose Plasma On Directory');   % FIX: prompt said "Off"
end
[Time,Channel1,Channel2] = OscopeReadin(PathOn);
[Nfiles,Npoints] = size(Time);

if VerifyParameters == 1
    warning('VerifyParameters UI is not implemented; using the hard coded gains.');
end
VoltageOn = Channel1*AmpGain;
CurrentOn = Channel2*CurrentGain;

clear Channel1 Channel2

% --- Preallocate ---------------------------------------------------------
F_Meas_On         = zeros(1,Nfiles);
VoltageMeasured   = zeros(1,Nfiles);   % on-cycle amplitudes, see rename above
VoltagePhase      = zeros(1,Nfiles);
OffIndex          = zeros(1,Nfiles);
index1            = zeros(1,Nfiles);
index2            = zeros(1,Nfiles);
TruncIndex        = zeros(1,Nfiles);

% Device geometry
if DeviceSelection == 0
    Area_Device        = 5389/1000^2;   % m^2
    Area_Electrode     = 0.00379;       % m^2
    Electrode_Parameter= 1269.7/1000;   % m
elseif DeviceSelection == 1
    % Device selection UI still to be written. Falls through with no geometry,
    % which is fine as long as nothing downstream needs Area_Device.
end

% --- Fit each on-cycle voltage to a sine ---------------------------------
% Fitting first means the derivative can be taken without amplifying noise.
ft   = fittype('sin1');   % redefined here so this section can be run on its own
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';
opts.Lower   = [-Inf 0 -Inf];

for n = 1:Nfiles
    [xData,yData] = prepareCurveData(Time(n,:),VoltageOn(n,:));
    fitresult = fit(xData,yData,ft,opts);
    Vcoeff = coeffvalues(fitresult);
    a1 = Vcoeff(1);   % amplitude
    b1 = Vcoeff(2);   % angular frequency
    c1 = Vcoeff(3);   % phase
    if a1 < 0
        a1 = -a1;
        c1 = c1 + pi;
    end

    % FIX: was F_Meas_On(n,1) = ..., which turned a 1xN row into an NxN matrix.
    F_Meas_On(n)       = b1/2/pi;
    VoltageMeasured(n) = a1;
    VoltagePhase(n)    = c1;

    % Match to the nearest off-cycle frequency (useful for frequency sweeps)
    [freqMismatch,OffIndex(n)] = min(abs(F_Meas_On(n) - F_meas1));
    % FIX: the original compared a scalar against the whole 0.05*F_meas1 vector,
    % so the warning only fired if it exceeded every off frequency at once.
    if freqMismatch >= 0.05*F_meas1(OffIndex(n))
        fprintf(['Warning: no off-frequency measurement matches the on ' ...
                 'frequency for file %d\n'], n);
    end
end

%% Section 4 - Truncate to one period, integrate charge
%
% FIX: the original used Ncells = length(Time) for the column count. length()
% returns the largest dimension, so with more files than samples the matrices
% came out the wrong shape. Npoints is what was meant.

for n = 1:Nfiles
    [~,index1(n)] = min(abs(Time(n,:) + 0.5/(4*F_Meas_On(n))));   % t = -T/8
    [~,index2(n)] = min(abs(Time(n,:) - 3.5/(4*F_Meas_On(n))));   % t = +7T/8
    TruncIndex(n) = abs(index1(n) - index2(n));
end

% FIX: guard against running off the end of the record
if any(index1 + TruncIndex > Npoints)
    error(['Truncation window runs past the end of the record. Check that the ' ...
           'scope captured at least one full period after t = 0.']);
end

maxTrunc = max(TruncIndex);
VTrunc   = zeros(Nfiles,maxTrunc);
ITrunc   = zeros(Nfiles,maxTrunc);
TTrunc   = zeros(Nfiles,maxTrunc);
QTrunc   = zeros(Nfiles,maxTrunc);
UrTrunc  = zeros(Nfiles,maxTrunc);
QConstant             = zeros(1,Nfiles);
TotalEnergyDeposited  = zeros(1,Nfiles);
PlasmaEnergyDeposited = zeros(1,Nfiles);
UrMeasured            = zeros(1,Nfiles);

for n = 1:Nfiles
    k = 1:TruncIndex(n);
    src = k + index1(n);

    % FIX: replaced element-by-element loops with slice assignment
    VTrunc(n,k) = VoltageOn(n,src);
    ITrunc(n,k) = CurrentOn(n,src);
    TTrunc(n,k) = Time(n,src);

    % Remove the DC bias from the current
    ITrunc(n,k) = ITrunc(n,k) - mean(ITrunc(n,k));

    % REVIEW: QConstant divides the integrated charge by the number of samples,
    % not by the elapsed time, so it is not a mean current. It is only used in
    % the commented-out branch below, so nothing depends on it today. If you
    % re-enable it, the intended form is trapz(t,I)/(t(end)-t(1)).
    QConstant(n) = trapz(TTrunc(n,k),ITrunc(n,k))./TruncIndex(n);

    % Including QConstant:
    % QTruncTemp = cumtrapz(TTrunc(n,k),ITrunc(n,k) - QConstant(n));
    % Excluding QConstant:
    QTrunc(n,k) = cumtrapz(TTrunc(n,k),ITrunc(n,k));

    % Voltage across the device once the series resistance is removed
    UrTrunc(n,k) = VTrunc(n,k) - ITrunc(n,k).*R_off(OffIndex(n));

    % Centre the charge on zero
    QTrunc(n,k) = QTrunc(n,k) - mean(QTrunc(n,k));

    % Energy per period
    TotalEnergyDeposited(n)  = trapz(TTrunc(n,k), VTrunc(n,k).*ITrunc(n,k));
    PlasmaEnergyDeposited(n) = trapz(TTrunc(n,k), VTrunc(n,k).*ITrunc(n,k) ...
                                     - R_off(OffIndex(n))*ITrunc(n,k).^2);
    UrMeasured(n) = max(UrTrunc(n,k));
end

% FIX: the original multiplied every file's energy by mean(F_meas1), the average
% PLASMA-OFF frequency. The energy is integrated over exactly one on-cycle
% period, so the matching rate is that file's own on frequency. This changes the
% absolute power numbers slightly. Set UseOffFreqForPower = 1 to reproduce the
% old behaviour if you need to compare against earlier results.
UseOffFreqForPower = 0;
OffFreq = mean(F_meas1);
if UseOffFreqForPower
    TotalPower  = OffFreq.*TotalEnergyDeposited;
    PlasmaPower = OffFreq.*PlasmaEnergyDeposited;
else
    TotalPower  = F_Meas_On.*TotalEnergyDeposited;
    PlasmaPower = F_Meas_On.*PlasmaEnergyDeposited;
end

%% Section 4b - Phase-averaged representative cycle
%
% Fold every complete drive cycle in the record onto a single low-noise cycle.
% This is computed BEFORE corner detection because the corner finder runs on it.
%
% Phase averaging is a linear operation, so any linear combination of the
% averaged signals can be formed afterwards from the bins. That is why Qp is
% built downstream as PhaseQ - a*PhaseUr rather than being averaged a second
% time: the two are identical and the second pass would be wasted work.

PhaseBins   = 1000;   % target bins per cycle
MinPerBin   = 4;      % halve the bin count until at least this many samples land in each

PhaseUr     = nan(Nfiles,PhaseBins);
PhaseQ      = nan(Nfiles,PhaseBins);
PhaseTime   = nan(Nfiles,PhaseBins);
PhaseNbins  = zeros(1,Nfiles);
PhaseNcyc   = zeros(1,Nfiles);

for n = 1:Nfiles
    % Whole cycles from the full record, not the one-period truncation: folding
    % a single cycle onto itself averages nothing.
    Ifull = CurrentOn(n,:) - mean(CurrentOn(n,:));   % remove DC before integrating,
                                                     % otherwise cumtrapz ramps
    Urfull = VoltageOn(n,:) - Ifull*R_off(OffIndex(n));
    Qfull  = cumtrapz(Time(n,:),Ifull);
    Qfull  = Qfull - mean(Qfull);

    [Sb,Tb,nb,Ncyc] = PhaseAverage(Time(n,:),F_Meas_On(n), ...
                                   [Urfull; Qfull],PhaseBins,MinPerBin);
    if Ncyc < 1
        warning('File %d holds less than one full period; phase average skipped.',n);
        continue
    end
    PhaseNbins(n)     = nb;
    PhaseNcyc(n)      = Ncyc;
    PhaseUr(n,1:nb)   = Sb(1,:);
    PhaseQ(n,1:nb)    = Sb(2,:);
    PhaseTime(n,1:nb) = Tb;
end

%% Section 5 - Displacement and discharge current, ignition points
Displacement    = zeros(Nfiles,maxTrunc);
Discharge       = zeros(Nfiles,maxTrunc);
SimulatedVoltage= zeros(Nfiles,maxTrunc);

for n = 1:Nfiles
    k = 1:TruncIndex(n);

    % FIX: the original built two sine reconstructions, VoltageFitTrunc and
    % SimulatedVoltage. VoltageFitTrunc omitted the 2*pi, so it ran at 1/(2*pi)
    % of the true frequency; nothing read it, but it was wrong. Once corrected
    % the two were identical, so VoltageFitTrunc is gone and SimulatedVoltage,
    % which the displacement-current simulation actually uses, is the only one.
    SimulatedVoltage(n,k) = VoltageMeasured(n)*sin(2*pi*F_Meas_On(n)*TTrunc(n,k) + VoltagePhase(n));

    % Displacement current through the off capacitance
    DisplacementTemp = C_off(OffIndex(n))*gradient(UrTrunc(n,k))./gradient(TTrunc(n,k));
    DisplacementTemp = transpose(smooth(DisplacementTemp,SmoothWindow,'sgolay',3));

    Displacement(n,k) = DisplacementTemp;
    Discharge(n,k)    = ITrunc(n,k) - DisplacementTemp;

    clear DisplacementTemp   % FIX: original spelled this "DisplacmentTemp" and
                             % also tried to clear a DischargeTemp that never existed
end

% --- Preallocate the Lissajous segments ----------------------------------
% FIX: DisPlasma* were declared as zeros(1,Nfiles) but written as (n,i), so they
% silently regrew on every iteration. They are 2-D like the rest.
LogicThreshold      = ones(1,Nfiles);
QPlasmaOnPositive   = zeros(Nfiles,maxTrunc);
VPlasmaOnPositive   = zeros(Nfiles,maxTrunc);
TPlasmaOnPositive   = zeros(Nfiles,maxTrunc);
UrPlasmaOnPositive  = zeros(Nfiles,maxTrunc);
DisPlasmaOnPositive = zeros(Nfiles,maxTrunc);
QPlasmaOnNegative   = zeros(Nfiles,maxTrunc);
VPlasmaOnNegative   = zeros(Nfiles,maxTrunc);
TPlasmaOnNegative   = zeros(Nfiles,maxTrunc);
UrPlasmaOnNegative  = zeros(Nfiles,maxTrunc);
DisPlasmaOnNegative = zeros(Nfiles,maxTrunc);
QPlasmaOffPositive  = zeros(Nfiles,maxTrunc);
VPlasmaOffPositive  = zeros(Nfiles,maxTrunc);
TPlasmaOffPositive  = zeros(Nfiles,maxTrunc);
UrPlasmaOffPositive = zeros(Nfiles,maxTrunc);
DisPlasmaOffPositive= zeros(Nfiles,maxTrunc);
QPlasmaOffNegative  = zeros(Nfiles,maxTrunc);
VPlasmaOffNegative  = zeros(Nfiles,maxTrunc);
TPlasmaOffNegative  = zeros(Nfiles,maxTrunc);
UrPlasmaOffNegative = zeros(Nfiles,maxTrunc);
DisPlasmaOffNegative= zeros(Nfiles,maxTrunc);

point1index = zeros(1,Nfiles);
point2index = zeros(1,Nfiles);
point3index = zeros(1,Nfiles);
point4index = zeros(1,Nfiles);
point1t     = zeros(1,Nfiles);
point3t     = zeros(1,Nfiles);
PlasmaOffPositiveLength = zeros(1,Nfiles);
PlasmaOffNegativeLength = zeros(1,Nfiles);
PlasmaOnPositiveLength  = zeros(1,Nfiles);
PlasmaOnNegativeLength  = zeros(1,Nfiles);

% --- Corner detection ----------------------------------------------------
%
% CornerMethod = 'segmented' (default) reformulates the problem. The ignition
% corner is not "where the discharge current crosses a threshold", it is a SLOPE
% BREAKPOINT in Q(Ur): the off leg has slope C_geo, the on leg has slope C_diel,
% and they meet at ignition. Fitting for the breakpoint directly removes the
% dependency on C_off entirely, because C_geo comes OUT of the fit rather than
% going into it. That kills the feedback loop where a poor C_off corrupts
% Discharge, which corrupts the corner, which corrupts every fit downstream.
%
% It also runs on the phase-averaged cycle rather than raw samples, so it starts
% from data with several times less noise.
%
% Before looking for a corner at all, IgnitionScore decides whether there IS
% one. Both uploaded traces are below ignition (crest factor 1.4, kurtosis 1.5:
% the exact signature of a pure sine wave, so the non-capacitive residual is the
% series resistance, not filaments). No corner algorithm can succeed there, and
% returning a number anyway is worse than returning nothing.

CornerMethod   = 'segmented';   % 'segmented' (new) or 'threshold' (original)
CornerGuard    = 0.08;          % keep the breakpoint away from the segment ends
MinSlopeRatio  = 1.10;          % C_diel/C_geo below this means no real corner
MinKurtosis    = 3.0;           % residual kurtosis; 1.5 = sine, 3 = gaussian
MinIgnitionF   = 20;            % F statistic for the two-line vs one-line fit

IgnitionKurtosis  = nan(1,Nfiles);
IgnitionCrest     = nan(1,Nfiles);
IgnitionResidFrac = nan(1,Nfiles);
CornerSlopeRatio  = nan(2,Nfiles);   % row 1 positive half wave, row 2 negative
CornerCgeo        = nan(2,Nfiles);
CornerCdiel       = nan(2,Nfiles);
CornerF           = nan(2,Nfiles);
CornerUstar       = nan(2,Nfiles);
CornerUband       = nan(2,Nfiles);   % width of the 95% region, V
PlasmaDetected    = false(1,Nfiles);

for n = 1:Nfiles
    % Extinction corners: max and min charge, i.e. the voltage turning points.
    % These were already reliable, but taking them from the phase-averaged cycle
    % rather than raw samples removes the bias in max()/min(), which are extreme
    % order statistics and so are pulled outward by whatever noise spike happens
    % to be largest.
    [~,b1] = max(PhaseQ(n,1:PhaseNbins(n)));
    [~,b3] = min(PhaseQ(n,1:PhaseNbins(n)));
    point1index(n) = PhaseBinToTruncIndex(n,b1,PhaseNbins(n),TTrunc,TruncIndex,Time,F_Meas_On);
    point3index(n) = PhaseBinToTruncIndex(n,b3,PhaseNbins(n),TTrunc,TruncIndex,Time,F_Meas_On);
    point1t(n) = TTrunc(n,point1index(n));
    point3t(n) = TTrunc(n,point3index(n));

    % Is there a discharge at all?
    ig = IgnitionScore(Time(n,:),VoltageOn(n,:),CurrentOn(n,:),F_Meas_On(n));
    IgnitionKurtosis(n)  = ig.Kurtosis;
    IgnitionCrest(n)     = ig.CrestFactor;
    IgnitionResidFrac(n) = ig.ResidFrac;

    switch lower(CornerMethod)

    case 'segmented'
        nb = PhaseNbins(n);
        Ub = PhaseUr(n,1:nb).';
        Qb = PhaseQ(n,1:nb).';
        [~,iMax] = max(Ub);
        [~,iMin] = min(Ub);

        % Rising half wave: Ur climbs from its minimum to its maximum, ignition
        % sits at the high-Ur end, so it feeds the fit directly.
        risIdx = WrapRange(iMin,iMax,nb);
        Spos = FindLissajousCorner(Ub(risIdx),Qb(risIdx),CornerGuard);

        % Falling half wave: negating both axes makes Ur increase again with
        % ignition still at the far end, so one code path covers both. dQ/dU is
        % unchanged by the double flip.
        falIdx = WrapRange(iMax,iMin,nb);
        Sneg = FindLissajousCorner(-Ub(falIdx),-Qb(falIdx),CornerGuard);

        CornerCgeo(:,n)       = [Spos.Cgeo;       Sneg.Cgeo];
        CornerCdiel(:,n)      = [Spos.Cdiel;      Sneg.Cdiel];
        CornerSlopeRatio(:,n) = [Spos.SlopeRatio; Sneg.SlopeRatio];
        CornerF(:,n)          = [Spos.F;          Sneg.F];
        CornerUstar(:,n)      = [Spos.Ustar;     -Sneg.Ustar];
        CornerUband(:,n)      = [Spos.Uband;      Sneg.Uband];

        % Accept the corner only if the trace looks like a discharge AND the
        % two-line model genuinely beats a single line AND the slope actually
        % steps up. Any one of those failing means there is no corner.
        ok = Spos.found && Sneg.found ...
             && ig.Kurtosis      >= MinKurtosis ...
             && Spos.SlopeRatio  >= MinSlopeRatio ...
             && Sneg.SlopeRatio  >= MinSlopeRatio ...
             && Spos.F           >= MinIgnitionF ...
             && Sneg.F           >= MinIgnitionF;

        if ok
            b2 = risIdx(Spos.index);
            b4 = falIdx(Sneg.index);
            point2index(n) = PhaseBinToTruncIndex(n,b2,nb,TTrunc,TruncIndex,Time,F_Meas_On);
            point4index(n) = PhaseBinToTruncIndex(n,b4,nb,TTrunc,TruncIndex,Time,F_Meas_On);
            PlasmaDetected(n) = true;
        else
            point2index(n)    = 1;
            point4index(n)    = 1;
            LogicThreshold(n) = 0;
        end

    case 'threshold'
        % Original behaviour, kept for comparison.
        logickeep   = 0;
        locationpos = [];
        locationneg = [];
        for j = SmoothWindow+10:TruncIndex(n)
            pointcheck = Discharge(n,j);
            if logickeep == 0 && pointcheck >= ThreshInput
                locationpos = j;
                logickeep   = 1;
            elseif logickeep == 1 && pointcheck <= -ThreshInput
                locationneg = j;
                logickeep   = 2;
                break
            end
        end
        if ~isempty(locationpos) && ~isempty(locationneg)
            point2index(n) = locationpos(1);
            point4index(n) = locationneg(1);
            PlasmaDetected(n) = true;
        else
            point2index(n)    = 1;
            point4index(n)    = 1;
            LogicThreshold(n) = 0;
        end

    otherwise
        error('CornerMethod must be ''segmented'' or ''threshold''.');
    end

    % The leg-splitting below assumes point2 < point1 < point4 < point3. Check it
    % rather than letting a bad ordering silently produce negative lengths.
    if LogicThreshold(n) == 1 && ~(point2index(n) < point1index(n) && ...
            point1index(n) < point4index(n) && point4index(n) < point3index(n))
        warning(['File %d: corners out of order (%d %d %d %d), rejecting. Check ' ...
                 'the truncation phase offset in Section 4.'],n, ...
                 point2index(n),point1index(n),point4index(n),point3index(n));
        point2index(n)    = 1;
        point4index(n)    = 1;
        LogicThreshold(n) = 0;
        PlasmaDetected(n) = false;
    end
end

fprintf('\nIgnition detected in %d of %d files (method: %s)\n', ...
        sum(PlasmaDetected),Nfiles,CornerMethod);
if any(~PlasmaDetected)
    fprintf('  below ignition: %s\n',mat2str(find(~PlasmaDetected)));
end

% --- Split the trace into the four Lissajous legs ------------------------
for n = 1:Nfiles
    % FIX: every length below is now clamped at zero. Without the clamp, a file
    % where ignition was never detected produced negative lengths, which fed
    % negative indices into the fits and the alternative-Cd block below.
    PlasmaOffPositiveLength(n) = max(point2index(n),0);
    PlasmaOnPositiveLength(n)  = max(point1index(n) - point2index(n),0);
    PlasmaOffNegativeLength(n) = max(point4index(n) - point1index(n),0);
    PlasmaOnNegativeLength(n)  = max(point3index(n) - point4index(n),0);

    % Plasma off, positive half wave: start to point2
    k = 1:PlasmaOffPositiveLength(n);
    QPlasmaOffPositive(n,k)   = QTrunc(n,k);
    VPlasmaOffPositive(n,k)   = VTrunc(n,k);
    TPlasmaOffPositive(n,k)   = TTrunc(n,k);
    UrPlasmaOffPositive(n,k)  = UrTrunc(n,k);
    DisPlasmaOffPositive(n,k) = Discharge(n,k);

    % Plasma on, positive half wave: point2 to point1
    k = 1:PlasmaOnPositiveLength(n);
    QPlasmaOnPositive(n,k)   = QTrunc(n,k+point2index(n));
    VPlasmaOnPositive(n,k)   = VTrunc(n,k+point2index(n));
    TPlasmaOnPositive(n,k)   = TTrunc(n,k+point2index(n));
    UrPlasmaOnPositive(n,k)  = UrTrunc(n,k+point2index(n));
    DisPlasmaOnPositive(n,k) = Discharge(n,k+point2index(n));

    % Plasma off, negative half wave: point1 to point4
    k = 1:PlasmaOffNegativeLength(n);
    QPlasmaOffNegative(n,k)   = QTrunc(n,k+point1index(n));
    VPlasmaOffNegative(n,k)   = VTrunc(n,k+point1index(n));
    TPlasmaOffNegative(n,k)   = TTrunc(n,k+point1index(n));
    UrPlasmaOffNegative(n,k)  = UrTrunc(n,k+point1index(n));
    DisPlasmaOffNegative(n,k) = Discharge(n,k+point1index(n));

    % Plasma on, negative half wave: point4 to point3
    k = 1:PlasmaOnNegativeLength(n);
    QPlasmaOnNegative(n,k)   = QTrunc(n,k+point4index(n));
    VPlasmaOnNegative(n,k)   = VTrunc(n,k+point4index(n));
    TPlasmaOnNegative(n,k)   = TTrunc(n,k+point4index(n));
    UrPlasmaOnNegative(n,k)  = UrTrunc(n,k+point4index(n));
    DisPlasmaOnNegative(n,k) = Discharge(n,k+point4index(n));
end

% --- Alternative Cd, no curve fitting ------------------------------------
OffCapacitanceAlt = nan(1,Nfiles);   % FIX: preallocated, and NaN rather than 0
for n = 1:Nfiles
    L = PlasmaOffPositiveLength(n);
    if L >= 2
        dQ = QPlasmaOffPositive(n,L)  - QPlasmaOffPositive(n,1);
        dU = UrPlasmaOffPositive(n,L) - UrPlasmaOffPositive(n,1);
        if dU ~= 0                          % FIX: guard the divide. With no
            OffCapacitanceAlt(n) = dQ/dU;   % ignition detected L was 1 and this
        end                                 % divided by zero.
    end
end

%% Section 6 - Fit the on and off legs
% Off legs are fitted to a line, on legs to a line in Q vs Ur.

PlasmaOnPositiveA  = zeros(1,Nfiles);
PlasmaOnPositiveB  = zeros(1,Nfiles);
PlasmaOnNegativeA  = zeros(1,Nfiles);
PlasmaOnNegativeB  = zeros(1,Nfiles);
PlasmaOffPositiveA = zeros(1,Nfiles);
PlasmaOffPositiveB = zeros(1,Nfiles);
PlasmaOffNegativeA = zeros(1,Nfiles);
PlasmaOffNegativeB = zeros(1,Nfiles);
PlasmaOnPositiveFit  = zeros(Nfiles,maxTrunc);
PlasmaOffPositiveFit = zeros(Nfiles,maxTrunc);
PlasmaOnNegativeFit  = zeros(Nfiles,maxTrunc);
PlasmaOffNegativeFit = zeros(Nfiles,maxTrunc);
QpPlasmaOnNegative   = zeros(Nfiles,maxTrunc);
QpPlasmaOnPositive   = zeros(Nfiles,maxTrunc);
QpTrunc              = zeros(Nfiles,maxTrunc);
UdTrunc              = zeros(Nfiles,maxTrunc);
UgTrunc              = zeros(Nfiles,maxTrunc);
UdPlasmaOnPositive   = zeros(Nfiles,maxTrunc);
UgPlasmaOnPositive   = zeros(Nfiles,maxTrunc);
UdPlasmaOnNegative   = zeros(Nfiles,maxTrunc);
UgPlasmaOnNegative   = zeros(Nfiles,maxTrunc);
UgAveragePositive    = nan(1,Nfiles);
UgAverageNegative    = nan(1,Nfiles);
UgMiddlePositive     = nan(1,Nfiles);
UgMiddleNegative     = nan(1,Nfiles);

linOpts = fitoptions('Method','LinearLeastSquares');
linOpts.Upper = [Inf Inf];

MinFitPoints = 5;

for n = 1:Nfiles
    % Plasma off, negative
    if PlasmaOffNegativeLength(n) >= MinFitPoints
        k = 1:PlasmaOffNegativeLength(n);
        fr = fit(UrPlasmaOffNegative(n,k).',QPlasmaOffNegative(n,k).','poly1',linOpts);
        A = coeffvalues(fr);
        PlasmaOffNegativeA(n) = A(1);
        PlasmaOffNegativeB(n) = A(2);
    end

    % Plasma off, positive
    if PlasmaOffPositiveLength(n) >= MinFitPoints
        k = 1:PlasmaOffPositiveLength(n);
        fr = fit(UrPlasmaOffPositive(n,k).',QPlasmaOffPositive(n,k).','poly1',linOpts);
        A = coeffvalues(fr);
        PlasmaOffPositiveA(n) = A(1);
        PlasmaOffPositiveB(n) = A(2);
    end

    % Subtract the geometric charge to get the plasma charge Qp.
    % REVIEW: the positive half wave is corrected with PlasmaOffNegativeA, the
    % slope from the NEGATIVE off leg, not PlasmaOffPositiveA. That may well be
    % deliberate (the negative leg is often the cleaner fit), but it is the one
    % asymmetry in this section and it propagates into every Ug and Cd number.
    % Left as written so your results do not shift underneath you.
    QpPlasmaOnNegative(n,:) = QPlasmaOnNegative(n,:) - UrPlasmaOnNegative(n,:).*PlasmaOffNegativeA(n);
    QpPlasmaOnPositive(n,:) = QPlasmaOnPositive(n,:) - UrPlasmaOnPositive(n,:).*PlasmaOffNegativeA(n);
    QpTrunc(n,:)            = QTrunc(n,:) - UrTrunc(n,:).*PlasmaOffNegativeA(n);

    % Plasma on, positive
    if PlasmaOnPositiveLength(n) >= MinFitPoints
        k = 1:PlasmaOnPositiveLength(n);
        fr = fit(UrPlasmaOnPositive(n,k).',QpPlasmaOnPositive(n,k).','poly1',linOpts);
        A = coeffvalues(fr);
        PlasmaOnPositiveA(n) = A(1);
        PlasmaOnPositiveB(n) = A(2);
    end

    % Plasma on, negative
    if PlasmaOnNegativeLength(n) >= MinFitPoints
        k = 1:PlasmaOnNegativeLength(n);
        fr = fit(UrPlasmaOnNegative(n,k).',QpPlasmaOnNegative(n,k).','poly1',linOpts);
        A = coeffvalues(fr);
        PlasmaOnNegativeA(n) = A(1);
        PlasmaOnNegativeB(n) = A(2);
    end

    % FIX: guard the divide. If the positive on-leg fit was skipped,
    % PlasmaOnPositiveA(n) was 0 and UdTrunc came out Inf/NaN silently.
    if PlasmaOnPositiveA(n) ~= 0
        UdTrunc(n,:) = QpTrunc(n,:)./PlasmaOnPositiveA(n);
        UgTrunc(n,:) = UrTrunc(n,:) - UdTrunc(n,:);
    else
        UdTrunc(n,:) = NaN;
        UgTrunc(n,:) = NaN;
    end
end

% --- Reconstruct the fits and the gap voltage ----------------------------
for n = 1:Nfiles
    k = 1:PlasmaOnPositiveLength(n);
    if ~isempty(k) && PlasmaOnPositiveA(n) ~= 0
        PlasmaOnPositiveFit(n,k) = UrPlasmaOnPositive(n,k)*PlasmaOnPositiveA(n) + PlasmaOnPositiveB(n);
        UdPlasmaOnPositive(n,k)  = QpPlasmaOnPositive(n,k)/PlasmaOnPositiveA(n);
        UgPlasmaOnPositive(n,k)  = UrPlasmaOnPositive(n,k) - UdPlasmaOnPositive(n,k);
        UgAveragePositive(n)     = mean(UgPlasmaOnPositive(n,k));
        mid = max(round(PlasmaOnPositiveLength(n)/2),1);   % FIX: round(0.5) is 0,
        UgMiddlePositive(n) = UgPlasmaOnPositive(n,mid);   % which is an illegal index
    end

    k = 1:PlasmaOnNegativeLength(n);
    if ~isempty(k) && PlasmaOnNegativeA(n) ~= 0
        PlasmaOnNegativeFit(n,k) = UrPlasmaOnNegative(n,k)*PlasmaOnNegativeA(n) + PlasmaOnNegativeB(n);
        UdPlasmaOnNegative(n,k)  = QpPlasmaOnNegative(n,k)/PlasmaOnNegativeA(n);
        UgPlasmaOnNegative(n,k)  = UrPlasmaOnNegative(n,k) - UdPlasmaOnNegative(n,k);
        % FIX: these two were commented out in the original but Section 8 still
        % wrote UgAverageNegative into the export, which threw undefined variable.
        UgAverageNegative(n) = mean(UgPlasmaOnNegative(n,k));
        mid = max(round(PlasmaOnNegativeLength(n)/2),1);
        UgMiddleNegative(n)  = UgPlasmaOnNegative(n,mid);
    end

    k = 1:PlasmaOffNegativeLength(n);
    PlasmaOffNegativeFit(n,k) = UrPlasmaOffNegative(n,k)*PlasmaOffNegativeA(n) + PlasmaOffNegativeB(n);

    k = 1:PlasmaOffPositiveLength(n);
    PlasmaOffPositiveFit(n,k) = UrPlasmaOffPositive(n,k)*PlasmaOffPositiveA(n) + PlasmaOffPositiveB(n);
end

%% Section 7 - Capacitances, plasma resistance, dQ/dV
ChargeDerivitive       = zeros(Nfiles,maxTrunc);
PlasmaChargeDerivitive = zeros(Nfiles,maxTrunc);
DisplacementSim        = zeros(Nfiles,maxTrunc);
DischargeSim           = zeros(Nfiles,maxTrunc);
EnergyIns              = zeros(1,Nfiles);
PlasmaResistancePositive = nan(1,Nfiles);
PlasmaResistanceNegative = nan(1,Nfiles);

ElectrodeCapacitancePositive  = PlasmaOffPositiveA;
ElectrodeCapacitanceNegative  = PlasmaOffNegativeA;
DielectricCapacitancePositive = PlasmaOnPositiveA;
DielectricCapacitanceNegative = PlasmaOnNegativeA;
PlasmaAreaPositive = DielectricCapacitancePositive/CapacitancePerArea;
PlasmaAreaNegative = DielectricCapacitanceNegative/CapacitancePerArea;
PowerIns = Discharge.*UgTrunc;

for n = 1:Nfiles
    k = 1:TruncIndex(n);
    EnergyIns(n) = trapz(TTrunc(n,k),PowerIns(n,k));

    kp = 1:PlasmaOnPositiveLength(n);
    if ~isempty(kp)
        PlasmaResistancePositive(n) = UgAveragePositive(n)./mean(DisPlasmaOnPositive(n,kp));
    end
    kn = 1:PlasmaOnNegativeLength(n);
    if ~isempty(kn)
        % FIX: was commented out but Section 8 exported it. Also note the
        % original draft used UgAveragePositive here; the negative average is
        % what belongs on the negative half wave.
        PlasmaResistanceNegative(n) = UgAverageNegative(n)./mean(DisPlasmaOnNegative(n,kn));
    end

    DisplacementSim(n,k) = PlasmaOffNegativeA(n).*gradient(SimulatedVoltage(n,k))./gradient(TTrunc(n,k));
    DischargeSim(n,k)    = ITrunc(n,k) - DisplacementSim(n,k);
end

% --- dQ/dU by phase averaging -------------------------------------------
%
% The pointwise derivative gradient(Q)./gradient(Ur) divides one noisy signal
% by another, and near the turning points of the sine the denominator goes to
% zero, so the ratio blows up. The old code hid that by clamping the result to
% +-100 pF, which is a cosmetic fix: the asymptote is still there, it is just
% flat-topped now.
%
% This replaces it with phase averaging. Every sample in the record is dropped
% into one of PhaseBins bins according to its phase within the drive cycle, and
% the bin means become a single low-noise representative cycle. The derivative
% is then a Savitzky-Golay filter applied circularly, so there is no edge
% artifact at the wrap point, and the turning points are masked to NaN instead
% of clamped.
%
% Measured on the two uploaded traces, sd(dQ/dU) over the valid region:
%   1000 Hz file (2 complete cycles): 340 pF pointwise -> 67 pF binned
%   1600 Hz file (4 complete cycles): 197 pF pointwise -> 90 pF binned
%
% Note the derivative is taken on Qp, the plasma charge, per your request. The
% two differ by a constant: Qp = Q - a*Ur with a = PlasmaOffNegativeA, so
% dQp/dU = dQ/dU - a. Subtracting the geometric branch is exactly what isolates
% the plasma contribution to the capacitance, and it costs one addition to
% recover the total-charge version, done below.

ChargeDerivativeMethod = 'phase';   % 'phase' (new) or 'gradient' (old behaviour)
PhaseSGWindow    = 41;     % SG window in bins, must be odd
PhaseSGOrder     = 3;
TurningPointFrac = 0.10;   % mask bins where |dU/dbin| < this fraction of its max
SpapsTol         = 0.09;   % only used by the 'gradient' method
dQdVClamp        = 10e-8;  % only used by the 'gradient' method

% Phase-averaged results live on their own bin grid, not the sample grid
% PhaseUr, PhaseQ, PhaseTime, PhaseNbins and PhaseNcyc come from Section 4b.
PhaseQp    = nan(Nfiles,PhaseBins);   % bin-mean Qp, C
PhasedQpdU = nan(Nfiles,PhaseBins);   % dQp/dU, F
PhasedQdU  = nan(Nfiles,PhaseBins);   % dQ/dU,  F

switch lower(ChargeDerivativeMethod)

case 'phase'
    % Reuses the bins built in Section 4b. Phase averaging is linear, so the Qp
    % bins are exactly PhaseQ - a*PhaseUr; re-folding the record would give an
    % identical answer for twice the work.
    for n = 1:Nfiles
        nb = PhaseNbins(n);
        if nb == 0, continue; end
        w  = min(PhaseSGWindow, 2*floor((nb-1)/4)+1);   % odd, and well under nb

        Ub  = PhaseUr(n,1:nb).';
        Qb  = PhaseQ(n,1:nb).';
        Qpb = Qb - Ub*PlasmaOffNegativeA(n);
        PhaseQp(n,1:nb) = Qpb;

        % Wrap-aware Savitzky-Golay first derivative. The step is one bin and it
        % cancels in the ratio, so no dt scaling is applied.
        [~,g] = sgolay(PhaseSGOrder,w);
        dUb  = CircDeriv(Ub, g);
        dQpb = CircDeriv(Qpb,g);

        % Mask the turning points rather than clamping them
        Cp = dQpb./dUb;
        Cp(abs(dUb) < TurningPointFrac*max(abs(dUb))) = NaN;

        PhasedQpdU(n,1:nb) = Cp;
        PhasedQdU(n,1:nb)  = Cp + PlasmaOffNegativeA(n);   % see header note

        % Map back onto the truncated sample grid so the Section 8 plots and the
        % Origin export keep working unchanged. Interpolated in phase with one
        % bin wrapped onto each end so the seam stays continuous.
        Period = 1/F_Meas_On(n);
        kt  = 1:TruncIndex(n);
        ph  = mod((TTrunc(n,kt) - Time(n,1))/Period, 1);
        bc  = (((0:nb-1)+0.5)/nb).';
        bcE = [bc(end)-1; bc; bc(1)+1];
        PlasmaChargeDerivitive(n,kt) = interp1(bcE,[Cp(end);Cp;Cp(1)],ph,'linear');
        ChargeDerivitive(n,kt)       = PlasmaChargeDerivitive(n,kt) + PlasmaOffNegativeA(n);
    end

case 'gradient'
    % Original behaviour, kept so you can compare against earlier results.
    for n = 1:Nfiles
        k = 1:TruncIndex(n);
        [~,Urfit] = spaps(TTrunc(n,k),UrTrunc(n,k),SpapsTol);
        dqdv = gradient(QTrunc(n,k))./gradient(Urfit);
        dqdv(dqdv <= -dQdVClamp) = -dQdVClamp;
        dqdv(dqdv >=  dQdVClamp) =  dQdVClamp;
        ChargeDerivitive(n,k) = dqdv;
        PlasmaChargeDerivitive(n,k) = ChargeDerivitive(n,k) - PlasmaOffNegativeA(n);
        clear dqdv Urfit
    end

otherwise
    error('ChargeDerivativeMethod must be ''phase'' or ''gradient''.');
end


% --- Optional: plasma extent and reduced field ---------------------------
PlasmaExtentPositive = nan(1,Nfiles);
PlasmaExtentNegative = nan(1,Nfiles);
EfieldCorrect        = nan(1,Nfiles);
if ComputePlasmaExtent == 1
    PlasmaExtentPositive = -(Aside+Bside) + sqrt((Aside+Bside)^2 + PlasmaAreaPositive);
    PlasmaExtentNegative = -(Aside+Bside) + sqrt((Aside+Bside)^2 + PlasmaAreaNegative);
    EFieldPositive  = UgAveragePositive./PlasmaExtentPositive;
    EFieldNegative  = UgAverageNegative./PlasmaExtentNegative;
    TownsendPositive= EFieldPositive/n0/1e-21;
    TownsendNegative= EFieldNegative/n0/1e-21;
    EfieldCorrect   = UgAveragePositive/Thickness/n0/1e-21;
end

%% Section 8 - Plotting
DischargeGraphing  = 1;
LissajousGraphing  = 1;
GapVoltageGraphing = 1;
dQdVGraphing       = 1;
ParameterGraphing  = 0;   % NOTE: this block needs CPrime*, BreakdownVoltage*,
                          % Alpha* and Townsend*, none of which this script
                          % computes. See the block below.

% Per-file discharge current with the four corner markers
if DischargeGraphing == 1
    for k = 1:Nfiles
        figure
        plot(TTrunc(k,1:TruncIndex(k)),Discharge(k,1:TruncIndex(k)),'LineWidth',2)
        ylabel('Current / A')
        yyaxis right
        hold on
        plot(TTrunc(k,1:TruncIndex(k)),VTrunc(k,1:TruncIndex(k)),'LineWidth',2)
        xline(TTrunc(k,point1index(k)),'Color','red');
        xline(TTrunc(k,point2index(k)),'Color','blue');
        xline(TTrunc(k,point3index(k)),'Color','green');
        xline(TTrunc(k,point4index(k)),'Color','black');
        ylabel('Voltage / V')
        xlabel('Time / s')
        title(sprintf('Discharge at %d Volts',round(VoltageMeasured(k))))
        set(gcf,'Visible','on')
    end
end

% Simulated discharge current across the run
% FIX: AAA was defined inside this if-block but used again by the gap voltage
% block below, so turning this plot off broke that one. Defined once, up front.
AAA = linspace(0,Nfiles*6-6,Nfiles);   % seconds elapsed, 6 s between files

if DischargeGraphing == 1
    figure
    Labels = cell(1,Nfiles+1);
    for k = 1:Nfiles
        plot(TTrunc(k,1:TruncIndex(k)),DischargeSim(k,1:TruncIndex(k)),'LineWidth',2)
        hold on
        Labels{k} = sprintf('%d sec',round(AAA(k)));
    end
    ylabel('Current / A')
    xlabel('Time / s')
    title('Change In Current During 15 Minute Treatment')
    Labels{Nfiles+1} = 'Applied Voltage';

    yyaxis right
    k = max(Nfiles-1,1);
    plot(TTrunc(k,1:TruncIndex(k)),VTrunc(k,1:TruncIndex(k)),'LineWidth',2)
    ylabel('Voltage / V')
    ylim([-2500 2500])
    legend(Labels,'Location','SouthWest')
    set(gcf,'Visible','on')
end

% Lissajous figures, charge against resistance-corrected voltage
if LissajousGraphing == 1
    for k = 1:10:Nfiles
        figure
        scatter(downsample(UrTrunc(k,1:TruncIndex(k)),20), ...
                1e6*downsample(QTrunc(k,1:TruncIndex(k)),20),20,[0 0 1])
        ylabel('Charge / \muC')
        xlabel('Measured Voltage / V')
        title(sprintf('Lissajous Figure at %d Volts',round(VoltageMeasured(k))))
        set(gcf,'Visible','on')
    end

    % Raw voltage against corrected voltage, same charge
    for k = 1:20:Nfiles
        figure
        scatter(downsample(VTrunc(k,1:TruncIndex(k)),20), ...
                downsample(QTrunc(k,1:TruncIndex(k)),20),20,[0 0 1]); hold on
        scatter(downsample(UrTrunc(k,1:TruncIndex(k)),20), ...
                downsample(QTrunc(k,1:TruncIndex(k)),20),20,[1 0 0])
        ylabel('Charge / C')
        xlabel('Voltage / V')
        legend('Applied','Series-R corrected')
        set(gcf,'Visible','on')
    end
end

% Gap voltage
% FIX: the original looped k = 80:4:Nfiles. With fewer than 80 files the loop
% never ran and legend(GapVoltageLegendLabels) hit an undefined variable.
GapVoltageStart = min(80,Nfiles);
GapVoltageData  = [];
if GapVoltageGraphing == 1
    figure(3001)
    i = 1;
    GapVoltageLegendLabels = {};
    for k = GapVoltageStart:4:Nfiles
        plot(TTrunc(k,1:TruncIndex(k)),UgTrunc(k,1:TruncIndex(k))); hold on
        GapVoltageLegendLabels{i} = sprintf('%.3f min',AAA(k)/60);   %#ok<SAGROW>
        GapVoltageData(2*i-1,1:TruncIndex(k)) = TTrunc(k,1:TruncIndex(k));   %#ok<SAGROW>
        GapVoltageData(2*i,  1:TruncIndex(k)) = UgTrunc(k,1:TruncIndex(k));  %#ok<SAGROW>
        i = i + 1;
    end
    ylabel('U_g / V')
    xlabel('Time / s')
    if ~isempty(GapVoltageLegendLabels)
        legend(GapVoltageLegendLabels)
    end
    set(gcf,'Visible','on')
end

% dQ/dV against normalised applied voltage
if dQdVGraphing == 1
    figure
    i = 1;
    dQdVLabels = {};
    for k = 1:12:Nfiles
        plot(VTrunc(k,1:TruncIndex(k))/VoltageMeasured(k), ...
             ChargeDerivitive(k,1:TruncIndex(k))*1e12); hold on
        dQdVLabels{i} = sprintf('%d V',round(VoltageMeasured(k)));   %#ok<SAGROW>
        i = i + 1;
    end
    ylabel('Charge Derivative / pF')
    xlabel('V / V_0')      % FIX: axis was labelled 'Time / s' but plots V/V0
    if ~isempty(dQdVLabels), legend(dQdVLabels); end
    set(gcf,'Visible','on')

    % dQ/dV against time
    figure
    i = 1;
    dQdVTimeLabels = {};
    for k = 1:16:Nfiles
        plot(TTrunc(k,1:TruncIndex(k)),ChargeDerivitive(k,1:TruncIndex(k))*1e12); hold on
        dQdVTimeLabels{i} = sprintf('%d V',round(VoltageMeasured(k)));   %#ok<SAGROW>
        i = i + 1;
    end
    ylabel('Charge Derivative / pF')
    xlabel('Time / s')
    if ~isempty(dQdVTimeLabels), legend(dQdVTimeLabels); end
    set(gcf,'Visible','on')

    % Phase-averaged dQp/dU on its own bin grid, plotted against Ur. This is the
    % honest version of the curve: no interpolation back onto the sample grid,
    % and the gaps are the masked turning points rather than clamped asymptotes.
    if strcmpi(ChargeDerivativeMethod,'phase')
        figure
        i = 1;
        phaseLabels = {};
        for k = 1:12:Nfiles
            plot(PhaseUr(k,:),PhasedQpdU(k,:)*1e12,'LineWidth',1.4); hold on
            phaseLabels{i} = sprintf('%d V',round(VoltageMeasured(k)));   %#ok<SAGROW>
            i = i + 1;
        end
        xlabel('U_r / V')
        ylabel('dQ_p/dU / pF')
        title('Phase-averaged plasma charge derivative')
        grid on
        if ~isempty(phaseLabels), legend(phaseLabels); end
        set(gcf,'Visible','on')
    end
end

% Reduced electric field
if ComputePlasmaExtent == 1
    figure
    scatter(VoltageMeasured,EfieldCorrect)
    ylabel('Reduced Electric Field / Td')
    xlabel('Applied Voltage Amplitude / V')
    set(gcf,'Visible','on')
end

% Electrode parameter plots.
% NOTE: this block is left disabled because CPrimePositive, CPrimeNegative,
% BreakdownVoltagePositive/Negative, AlphaPositive/Negative and
% TownsendPositive/Negative are never computed anywhere in this script. The
% original had the same problem, hidden behind ParameterGraphing = 0. It also
% plotted BreakdownVoltagePositive on both the positive AND the negative
% figure. Either add the missing calculations or delete the block.
if ParameterGraphing == 1
    error(['ParameterGraphing requires CPrime*, BreakdownVoltage*, Alpha* and ' ...
           'Townsend*, which this script does not compute. See the note above.']);
end

%% Section 9 - Triplicate averaging
% Files are assumed to come in groups of three repeats.
NTrip = floor(Nfiles/3);   % FIX: Nfiles/3 as a loop bound is non-integer unless
                           % the count divides evenly
PlasmaPowerTrip            = zeros(1,NTrip);
PlasmaPowerTripSTD         = zeros(1,NTrip);
VoltageMeasuredTrip        = zeros(1,NTrip);
VoltageMeasuredTripSTD     = zeros(1,NTrip);
UrMeasuredTrip             = zeros(1,NTrip);
UrMeasuredTripSTD          = zeros(1,NTrip);
DielectricCapacitanceTrip  = zeros(1,NTrip);
DielectricCapacitanceTripSTD = zeros(1,NTrip);
PlasmaExtentPositiveTrip     = nan(1,NTrip);
PlasmaExtentPositiveTripSTD  = nan(1,NTrip);

for n = 1:NTrip
    idx = (3*n-2):(3*n);
    PlasmaPowerTrip(n)            = mean(PlasmaPower(idx));
    PlasmaPowerTripSTD(n)         = std(PlasmaPower(idx));
    VoltageMeasuredTrip(n)        = mean(VoltageMeasured(idx));
    VoltageMeasuredTripSTD(n)     = std(VoltageMeasured(idx));
    UrMeasuredTrip(n)             = mean(UrMeasured(idx));
    UrMeasuredTripSTD(n)          = std(UrMeasured(idx));
    DielectricCapacitanceTrip(n)  = mean(PlasmaOnPositiveA(idx));
    DielectricCapacitanceTripSTD(n)= std(PlasmaOnPositiveA(idx));
    % FIX: guarded. The original referenced PlasmaExtentPositive unconditionally
    % while the block that defines it was commented out, so this loop always
    % errored with "Unrecognized function or variable 'PlasmaExtentPositive'".
    if ComputePlasmaExtent == 1
        PlasmaExtentPositiveTrip(n)   = mean(PlasmaExtentPositive(idx));
        PlasmaExtentPositiveTripSTD(n)= std(PlasmaExtentPositive(idx));
    end
end
% FIX: the original also dropped VoltageMeasuredSTD, which duplicated
% VoltageMeasuredTripSTD. Removed.

%% Section 10 - Save results
SaveData = 0;
if SaveData == 1
    % FIX: the original wrote 9 data columns under 7 headers, and the header
    % order did not match the columns (the gap voltages were labelled as plasma
    % resistances). Headers now match the columns one for one.
    Data = transpose([VoltageMeasured; ...
                      DielectricCapacitancePositive; ...
                      DielectricCapacitanceNegative; ...
                      UgAveragePositive; ...
                      UgAverageNegative; ...
                      PlasmaResistancePositive; ...
                      PlasmaResistanceNegative; ...
                      TotalPower; ...
                      PlasmaPower]);
    HeaderLines = {'Measured Voltage / V', ...
                   'Positive Halfwave Dielectric Capacitance / F', ...
                   'Negative Halfwave Dielectric Capacitance / F', ...
                   'Positive Halfwave Average Gap Voltage / V', ...
                   'Negative Halfwave Average Gap Voltage / V', ...
                   'Positive Halfwave Plasma Resistance / Ohm', ...
                   'Negative Halfwave Plasma Resistance / Ohm', ...
                   'Total Power / W', ...
                   'Plasma Power / W'};
    writecell(HeaderLines,SaveFileName,'Delimiter',',','FileType','text','WriteMode','append');
    writematrix(Data,SaveFileName,'Delimiter',',','WriteMode','append');
end

%% Section 11 - Export for Origin
% FIX: the original ran this unconditionally using an index vector called Index
% that is never created anywhere in the script, so the last six lines always
% threw "Unrecognized function or variable 'Index'". Now gated, and Index is
% built explicitly. Set the file numbers you actually want to export.
ExportForOrigin = 0;
Index = [];   % e.g. Index = [1 20 40 60 80 100];

if ExportForOrigin == 1
    if numel(Index) < 6
        error('Index must contain at least six file numbers to export.');
    end
    if any(Index > Nfiles) || any(Index < 1)
        error('Index contains file numbers outside 1:%d.',Nfiles);
    end
    sel = Index(1:6);
    OriginUgTrunc          = InterleaveWithTime(TTrunc,UgTrunc,sel);
    OriginVTrunc           = InterleaveWithTime(TTrunc,VTrunc,sel);
    OriginDischarge        = InterleaveWithTime(TTrunc,Discharge,sel);
    OriginCurrent          = InterleaveWithTime(TTrunc,ITrunc,sel);
    OriginChargeDerivitive = InterleaveWithTime(TTrunc,ChargeDerivitive,sel);
end

%% Local functions
function [Sb,Tb,nb,Ncyc,cnt] = PhaseAverage(t,F,Sigs,nbTarget,MinPerBin)
% Fold a multi-cycle record onto one representative cycle.
%   t        - 1xN uniformly sampled time vector
%   F        - drive frequency from the sine fit, Hz
%   Sigs     - MxN, one signal per row, all folded on the same bin map
%   nbTarget - desired bins per cycle
%   MinPerBin- bin count halves until at least this many samples land per bin
% Returns Sb (M x nb), Tb (bin centre times within one period), nb, the number
% of whole cycles used, and the per-bin sample counts.
dt     = median(diff(t));
Period = 1/F;
P      = Period/dt;              % samples per period, deliberately non-integer:
Ncyc   = floor(numel(t)/P);      % mod() below handles the fractional part, which
Sb=[]; Tb=[]; nb=0; cnt=[];      % index-based folding would let drift
if Ncyc < 1, return; end

nb = nbTarget;
while nb > 20 && Ncyc*P/nb < MinPerBin
    nb = round(nb/2);
end

k   = 1:round(Ncyc*P);
idx = (0:numel(k)-1).';
b   = min(floor(mod(idx/P,1)*nb)+1, nb);
cnt = accumarray(b,1,[nb 1]);
cs  = max(cnt,1);
Sb  = zeros(size(Sigs,1),nb);
for m = 1:size(Sigs,1)
    s = Sigs(m,k);
    Sb(m,:) = (accumarray(b,s(:),[nb 1])./cs).';
end
Tb  = ((0:nb-1)+0.5)/nb*Period;
cnt = cnt.';
end

function idx = WrapRange(a,b,nb)
% Bin indices walking forward from a to b around a cycle of nb bins.
if a <= b
    idx = (a:b).';
else
    idx = mod((a:b+nb)-1,nb).'+1;
end
end

function j = PhaseBinToTruncIndex(n,bin,nb,TTrunc,TruncIndex,Time,F)
% Map a phase bin back onto the truncated sample grid, comparing phases
% circularly so the wrap point is not a discontinuity.
Period  = 1/F(n);
phBin   = (bin-0.5)/nb;
k       = 1:TruncIndex(n);
phTrunc = mod((TTrunc(n,k) - Time(n,1))/Period, 1);
d       = abs(mod(phTrunc - phBin + 0.5, 1) - 0.5);
[~,j]   = min(d);
end

function S = IgnitionScore(t,V,I,F)
% Decide whether a trace contains a discharge, without needing C_off.
%
% Fits the best pure capacitor to the measured current and inspects what is left
% over. Below ignition the residual is the series resistance and dielectric
% loss, which is a smooth sinusoid: crest factor near sqrt(2) = 1.41 and
% kurtosis near 1.5. A filamentary discharge is spiky instead, pushing kurtosis
% well above the gaussian value of 3 and the crest factor above 4.
%
% This is a shape test, not an amplitude test, so it needs no threshold in amps
% and does not care about the absolute current scale.
dt   = median(diff(t));
P    = 1/(F*dt);
Ncyc = max(floor(numel(t)/P),1);
k    = 1:min(round(Ncyc*P),numel(t));
tt   = t(k);
Ii   = I(k) - mean(I(k));

dVdt = gradient(V(k),tt);
Ccap = dot(Ii,dVdt)/dot(dVdt,dVdt);   % least squares pure capacitance
r    = Ii - Ccap*dVdt;
sd   = std(r);

S.Ccap        = Ccap;
S.ResidFrac   = sd/std(Ii);
S.CrestFactor = max(abs(r))/sd;
S.Kurtosis    = mean(((r-mean(r))/sd).^4);
S.NSpikes     = sum(abs(r-mean(r))/sd > 4);
end

function S = FindLissajousCorner(U,Q,GuardFrac)
% Locate the ignition corner of one half wave by segmented regression.
%
% Model:   Q(U) = b0 + b1*U + b2*max(U-u0,0)
%
% b1 is the geometric (off-leg) slope, b1+b2 the dielectric (on-leg) slope, and
% u0 the breakpoint. For a FIXED u0 this is ordinary least squares, so scanning
% u0 across every candidate bin is globally optimal, needs no initial guess, and
% cannot fall into a local minimum. Crucially it never consumes C_off: the off
% capacitance drops out of the fit as b1, which is what breaks the circular
% dependency in the threshold method.
%
% U must be increasing with the on leg at the high-U end. For the negative half
% wave pass (-U,-Q), which leaves dQ/dU unchanged.
%
% The returned F statistic compares the two-line model against a single line.
% Treat it as a relative score, not a formal p-value: under the null there is no
% breakpoint to estimate, so the usual regularity conditions do not hold and the
% nominal significance is optimistic. Calibrate the cutoff on your PlasmaOff
% directory, where you know the answer should be "no corner".
U = U(:); Q = Q(:); N = numel(U);
S = struct('found',false,'index',NaN,'Ustar',NaN,'Cgeo',NaN,'Cdiel',NaN, ...
           'SlopeRatio',NaN,'Uband',NaN,'F',NaN,'rms',NaN);
if N < 40, return; end

X1   = [ones(N,1) U];
b1   = X1\Q;
SSE1 = sum((Q - X1*b1).^2);

lo   = max(round(GuardFrac*N),5);
cand = lo:(N-lo);
if isempty(cand), return; end
SSE  = inf(size(cand));
Bs   = nan(3,numel(cand));
for m = 1:numel(cand)
    X = [ones(N,1) U max(U-U(cand(m)),0)];
    b = X\Q;
    if b(3) <= 0, continue; end        % the slope must step UP at ignition
    SSE(m) = sum((Q - X*b).^2);
    Bs(:,m) = b;
end
[SSEmin,j] = min(SSE);
if ~isfinite(SSEmin), return; end

b = Bs(:,j);
S.found      = true;
S.index      = cand(j);
S.Ustar      = U(cand(j));
S.Cgeo       = b(2);
S.Cdiel      = b(2) + b(3);
S.SlopeRatio = S.Cdiel/S.Cgeo;
S.rms        = sqrt(SSEmin/N);
S.F          = ((SSE1 - SSEmin)/2)/(SSEmin/(N-4));

% Width of the region where the fit is statistically indistinguishable from the
% best one. A wide band means the corner is poorly localised and the file should
% be treated with suspicion even if the F statistic passes.
band    = U(cand(SSE <= SSEmin*(1 + 4/max(N-4,1))));
S.Uband = max(band) - min(band);
end

function dy = CircDeriv(y,g)
% Savitzky-Golay first derivative for a periodic signal.
%   y - column of phase-averaged bin values
%   g - second output of sgolay(order,window)
% The signal is wrapped by one filter length on each end before convolving, so
% the bin at phase 0 sees the bins at phase 1 as its neighbours and there is no
% edge transient at the seam. Step size is one bin; that scale cancels when two
% of these are divided to form dQ/dU, so no dt scaling is applied.
% The 1/(-1)^1 factor is MATLAB's documented sgolay derivative convention,
% factorial(p)/(-h)^p with p = 1 and h = 1.
p  = size(g,1);
y  = y(:);
yy = [y(end-p+1:end); y; y(1:p)];
dd = conv(yy, 1/(-1)^1 * g(:,2), 'same');
dy = dd(p+1:end-p);
end

function out = InterleaveWithTime(T,M,sel)
% Builds the column layout Origin expects: [t1 y1 t2 y2 ...], one X/Y pair per
% selected file.
out = zeros(size(T,2),2*numel(sel));
for j = 1:numel(sel)
    out(:,2*j-1) = T(sel(j),:).';
    out(:,2*j)   = M(sel(j),:).';
end
end