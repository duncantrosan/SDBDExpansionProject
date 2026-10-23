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
OptionalInputs   = 1;
VerifyParameters = 0;
SingleFrequency  = 1;
DeviceSelection  = 1;
BolsigResults    = '';

% --- Numerical parameters ------------------------------------------------
SmoothWindow      = 41;
ThreshInput       = 0.00085;   % A (0.35 mA)
Perimeter         = 0.2266;    % m
CapacitancePerArea= 1.245e-7;  % F/m^2
Thickness         = 90e-6;     % m

% --- Probe / amplifier calibration ---------------------------------------
AmpGain       = 500;      % V/V
CurrentGain   = 6/1000;   % A/V (1 V = 6 mA)

% --- Loss tangent model parameters ---------------------------------------
% These come from the Hioki impedance analyser measurements on the device
% alone (no amplifier). If you have Hioki data, run HiokiImpedanceSweep.m
% first to get these numbers for your specific device.
%
% If you do not have Hioki data, leave UseHiokiTanDelta = 0 and the code
% will extract everything from the scope data alone (the old method). The
% results will be correct but R_esr and tan(d) will be coupled and harder
% to separate.
UseHiokiTanDelta = 0;          % 1 = use the known tan(d) below
TanDelta_Hioki   = 0.019;      % from Hioki: device loss tangent (constant with f)
                               % Kapton alone is 0.002; the excess is the adhesive.

% --- Optional geometry for plasma extent ---------------------------------
ComputePlasmaExtent = 1;
Aside = 0.0896;    % m
Bside = 0.0237;    % m
n0    = 2.686e25;   % m^-3

% --- Paths ---------------------------------------------------------------
PathOff = 'E:\MATLAB\Data\SDBDModelingProject\VoltageandCurrent\2025\1600 Hz\bighex\bighex\freqsweep\PlasmaOff';
PathOn  = 'E:\MATLAB\Data\SDBDModelingProject\VoltageandCurrent\2025\1600 Hz\bighex\bighex\vsweep\PlasmaOn';
DefaultBrowseDir = pwd;
RunName      = 'Hex Old';
SaveFileName = fullfile(pwd,'HexLattice.csv');

% --- Read the plasma-off traces ------------------------------------------
if isempty(PathOff)
    PathOff = uigetdir(DefaultBrowseDir,'Choose Plasma Off Directory');
end
[Time,Channel1,Channel2] = OscopeReadin(PathOff);
[Nfiles,Npoints] = size(Time);

if VerifyParameters == 1
    warning('VerifyParameters UI is not implemented; using hard coded gains.');
end
Voltage = Channel1*AmpGain;
Current = Channel2*CurrentGain;
clear Channel1 Channel2

%% Section 2 - Off measurements: R_esr, C, and tan(d)
%
% The device is modelled as:
%
%     R_esr ─── [ C with loss tangent tan(d) ]
%
% where R_esr is the EXTERNAL series resistance (amplifier output impedance,
% cable, connections) and the capacitor carries a constant loss tangent from
% the dielectric + adhesive stack.
%
% The total series-equivalent impedance at angular frequency w is:
%
%     Z = R_esr + tan(d)/(w*C) - j/(w*C)
%
% The old code reported a single "R_off" which was Re(Z) = R_esr + tan(d)/(wC).
% That number changes with frequency even though nothing physical is changing,
% because it mixes a constant (R_esr) with a 1/f term (dielectric loss).
%
% This version separates them. R_esr is what you subtract from V to get Ur
% (the voltage across the device terminals). The dielectric loss stays inside
% the device where it belongs.
%
% KEY CONSEQUENCE FOR THE REST OF THE CODE:
%   Old:  Ur = V - I * R_off           <- over-corrects by tan(d)/(wC)
%   New:  Ur = V - I * R_esr           <- correct
%
% At 1600 Hz with C = 480 pF and tan(d) = 0.019:
%   Old R_off  ~ 17,200 ohm  (R_esr + 4,100 ohm of dielectric loss)
%   New R_esr  ~ 13,200 ohm
%   Difference in Ur ~ 20 V on a 1029 V signal (2%)
%
% That 2% matters for Lissajous analysis because it tilts the parallelogram.

% --- Preallocate ---------------------------------------------------------
F_meas1            = zeros(1,Nfiles);
VoltageMeasuredOff = zeros(1,Nfiles);
I_meas             = zeros(1,Nfiles);
phi1               = zeros(1,Nfiles);
phi2               = zeros(1,Nfiles);
Z_off              = zeros(1,Nfiles);   % complex impedance per file
R_total            = zeros(1,Nfiles);   % Re(Z), the old "R_off"
C_off              = zeros(1,Nfiles);   % series capacitance
TanDelta_off       = zeros(1,Nfiles);   % measured tan(d) per file

ft   = fittype('sin1');
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';
opts.Lower   = [-Inf 0 -Inf];

for n = 1:Nfiles
    % --- Fit voltage to sine ---------------------------------------------
    [xData,yData] = prepareCurveData(Time(n,:),Voltage(n,:));
    fitresult = fit(xData,yData,ft,opts);
    Vcoeff = coeffvalues(fitresult);
    a1 = Vcoeff(1);  b1 = Vcoeff(2);  c1 = Vcoeff(3);
    if a1 < 0, a1 = -a1; c1 = c1 + pi; end
    F_meas1(n)            = b1/2/pi;
    VoltageMeasuredOff(n) = a1;
    phi1(n)               = c1;

    % --- Fit current to sine ---------------------------------------------
    [xData,yData] = prepareCurveData(Time(n,:),Current(n,:));
    fitresult = fit(xData,yData,ft,opts);
    Icoeff = coeffvalues(fitresult);
    a2 = Icoeff(1);  b2 = Icoeff(2);  c2 = Icoeff(3);
    if a2 < 0, a2 = -a2; c2 = c2 + pi; end
    I_meas(n) = a2;
    phi2(n)   = c2;

    % --- Complex impedance from phasor division --------------------------
    w  = F_meas1(n)*2*pi;
    vi = VoltageMeasuredOff(n)*exp(1i*(pi/2 - (w*Time(n,:) + phi1(n))));
    ii = I_meas(n)*exp(1i*(pi/2 - (w*Time(n,:) + phi2(n))));
    Z_off(n) = mean(vi./ii);

    R_total(n)     = real(Z_off(n));           % Re(Z) = R_esr + td/(wC)
    C_off(n)       = 1/(imag(Z_off(n))*w);    % series capacitance (sign convention)
    TanDelta_off(n)= R_total(n)*w*C_off(n);   % total tan(d) of the network
end

% --- Separate R_esr from dielectric loss ---------------------------------
%
% Method depends on what data you have:
%
% (A) If you have Hioki tan(d) from the device alone, subtract it:
%       R_esr = Re(Z) - tan(d)_Hioki / (w*C)
%     This is the most accurate because tan(d) is measured independently.
%
% (B) If you only have the scope sweep at multiple frequencies, fit:
%       Re(Z) = R_esr + tan(d)/(w*C)
%     Two unknowns, linear in 1/w, solved by least squares.
%
% (C) If you only have one frequency, R_esr and tan(d) are degenerate.
%     Report R_total and note the ambiguity.

if UseHiokiTanDelta == 1
    % Method A: known tan(d)
    TanDelta_device = TanDelta_Hioki;
    R_esr = R_total - TanDelta_device./(2*pi*F_meas1.*C_off);
    R_esr_mean = mean(R_esr);
    R_esr_err  = std(R_esr);
    fprintf('\nUsing Hioki tan(d) = %.4f\n',TanDelta_device);

elseif numel(unique(round(F_meas1))) >= 3
    % Method B: multi-frequency fit
    % Re(Z) = R_esr + (td/C) * 1/(2*pi*f)
    % Linear regression: y = R_total, x = 1/(2*pi*f*C)
    X_inv = 1./(2*pi*F_meas1.*C_off);          % 1/(wC) for each file
    A = [ones(Nfiles,1), X_inv(:)];
    beta = A \ R_total(:);
    R_esr_mean     = beta(1);
    TanDelta_device= beta(2);
    R_esr = R_esr_mean * ones(1,Nfiles);        % constant by definition
    R_esr_err = sqrt(mean((R_total(:) - A*beta).^2));  % residual as uncertainty
    fprintf('\nFitted from frequency sweep:\n');

else
    % Method C: single frequency, cannot separate
    warning(['Only one frequency in the off data. Cannot separate R_esr from ' ...
             'tan(d). Using total Re(Z) as R_off. Consider taking off ' ...
             'measurements at multiple frequencies.']);
    TanDelta_device = 0;
    R_esr = R_total;
    R_esr_mean = mean(R_esr);
    R_esr_err  = std(R_esr);
end

% --- Report results ------------------------------------------------------
C_off_mean = mean(C_off);
C_off_err  = std(C_off);
TanDelta_mean = mean(TanDelta_off);
TanDelta_err  = std(TanDelta_off);

% Equivalent parallel resistance at the mean frequency
f_mean = mean(F_meas1);
w_mean = 2*pi*f_mean;
if TanDelta_device > 0
    Rp_equiv = 1/(w_mean*C_off_mean*TanDelta_device);
else
    Rp_equiv = Inf;
end

fprintf('  R_esr (external series)    = %.1f +- %.1f Ohm  (%.2f kOhm)\n', ...
        R_esr_mean, R_esr_err, R_esr_mean/1e3);
fprintf('  C_off (device capacitance) = %.2f +- %.2f pF\n', ...
        C_off_mean*1e12, C_off_err*1e12);
fprintf('  tan(d) (device loss)       = %.5f\n', TanDelta_device);
fprintf('  Rp @ %.0f Hz (equivalent)  = %.2f MOhm\n', f_mean, Rp_equiv/1e6);
fprintf('  R_total = Re(Z) @ %.0f Hz  = %.1f Ohm  (old "R_off")\n', ...
        f_mean, mean(R_total));
fprintf('  Dielectric loss in R_total = %.1f Ohm  (%.1f%% of R_total)\n', ...
        TanDelta_device/(w_mean*C_off_mean), ...
        100*TanDelta_device/(w_mean*C_off_mean)/mean(R_total));
fprintf('\n');
fprintf('  Kapton alone: tan(d) = 0.002 (DuPont datasheet)\n');
fprintf('  Your device:  tan(d) = %.3f  -> %.0fx higher\n', ...
        TanDelta_device, TanDelta_device/0.002);
fprintf('  The excess is the adhesive layer.\n');

% --- Diagnostic plots ----------------------------------------------------
% R_total vs frequency (should show the 1/f trend from tan(d))
figure
scatter(F_meas1, R_total, 40, 'filled'); hold on
if TanDelta_device > 0
    fplot = linspace(min(F_meas1)*0.9, max(F_meas1)*1.1, 200);
    Rfit  = R_esr_mean + TanDelta_device./(2*pi*fplot*C_off_mean);
    plot(fplot, Rfit, 'r-', 'LineWidth', 1.5)
    yline(R_esr_mean, '--k', 'R_{esr}', 'LineWidth', 1.2)
    legend('Measured Re(Z)', ...
           sprintf('R_{esr} + tan(\\delta)/(\\omegaC), tan(\\delta) = %.4f', TanDelta_device), ...
           sprintf('R_{esr} = %.0f \\Omega', R_esr_mean), ...
           'Location', 'NorthEast')
end
xlabel('Frequency / Hz')
ylabel('Re(Z) / \Omega')
title('Series Resistance Decomposition')
grid on; set(gcf,'Visible','on')

% Capacitance vs trial
figure
bar(C_off*1e12)
xlabel('Trial Number'); ylabel('Capacitance / pF')
title('Off Capacitance'); set(gcf,'Visible','on')

% tan(d) vs frequency
figure
scatter(F_meas1, TanDelta_off, 40, 'filled')
if TanDelta_device > 0
    hold on; yline(TanDelta_device, '--r', 'Device tan(\delta)', 'LineWidth', 1.2)
end
xlabel('Frequency / Hz'); ylabel('tan \delta')
title('Total Network Loss Tangent vs Frequency')
grid on; set(gcf,'Visible','on')

% --- Store results for downstream use ------------------------------------
% R_esr is what the rest of the code should use for the Ur correction.
% CapacitanceAverage is used by the on-cycle analysis for Lissajous fitting.
ResistanceAverage  = R_esr_mean;   % USE THIS for Ur = V - I*R_esr
CapacitanceAverage = C_off_mean;
LossTangent        = TanDelta_device;
ZZ = Z_off;

% For backward compatibility, keep R_off as a vector but set it to R_esr
% so that every downstream "R_off(OffIndex(n))" picks up the correct value.
R_off = R_esr;

fprintf('NOTE: R_off is now R_esr (%.0f Ohm), not Re(Z) (%.0f Ohm).\n', ...
        R_esr_mean, mean(R_total));
fprintf('      Ur correction is smaller by %.0f Ohm (%.1f%% of V at 1 kHz).\n', ...
        mean(R_total)-R_esr_mean, ...
        100*(mean(R_total)-R_esr_mean)*w_mean*C_off_mean*VoltageMeasuredOff(end)/VoltageMeasuredOff(end));
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
IgnitionVoltagePositive = zeros(1,Nfiles);
IgnitionVoltageNegative = zeros(1,Nfiles);
for n = 1:Nfiles
    % The Lissajous figure is fitted as a parallelogram whose corners give the
    % plasma parameters. Assuming the scope settings the code expects:
    %   point1 = end of the plasma-off positive half wave (max charge)
    %   point2 = positive ignition
    %   point3 = end of the plasma-off negative half wave (min charge)
    %   point4 = negative ignition
    [~,point1index(n)] = max(QTrunc(n,1:TruncIndex(n)));   % FIX: was searching the
    [~,point3index(n)] = min(QTrunc(n,1:TruncIndex(n)));   % zero-padded tail too
    point1t(n) = TTrunc(n,point1index(n));                 % FIX: was indexing Time,
    point3t(n) = TTrunc(n,point3index(n));                 % not the truncated record

    % FIX: the original read point1x/point3x out of VoltageDrop, a matrix that
    % was allocated and never written, so both were always exactly zero. The
    % voltage at those corners is UrTrunc. Neither value is used downstream.
    point1x = UrTrunc(n,point1index(n));
    point3x = UrTrunc(n,point3index(n));

    % Threshold crossing on the discharge current marks ignition
    logickeep   = 0;
    locationpos = [];
    locationneg = [];
    % Start past the smoothing window so filter artifacts cannot trip the logic
    for j = SmoothWindow+10:TruncIndex(n)
        pointcheck = Discharge(n,j);
        if logickeep == 0 && pointcheck >= ThreshInput
            locationpos = j;
            logickeep   = 1;
        elseif logickeep == 1 && pointcheck <= -ThreshInput   % FIX: was a second
            locationneg = j;                                  % if, which could run
            logickeep   = 2;                                  % in the same iteration
            break                                             % FIX: stop once both found
        end
    end

    if ~isempty(locationpos) && ~isempty(locationneg)
        point2index(n) = locationpos(1);
        point4index(n) = locationneg(1);
    else
        point2index(n)    = 1;
        point4index(n)    = 1;
        LogicThreshold(n) = 0;
    end
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
    IgnitionVoltagePositive(n) = UrTrunc(n,point2index(n));
    IgnitionVoltageNegative(n) = UrTrunc(n,point4index(n));
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

    %DisplacementSim(n,k) = PlasmaOffNegativeA(n).*gradient(SimulatedVoltage(n,k))./gradient(TTrunc(n,k)) + C_off .* LossTangent .* 2*pi.*F_Meas_On(n) .* SimulatedVoltage(n,k) ;
    %DischargeSim(n,k)    = ITrunc(n,k) - DisplacementSim(n,k);
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
PhaseBins        = 1000;   % bins per drive cycle
PhaseSGWindow    = 41;     % SG window in bins, must be odd
PhaseSGOrder     = 3;
TurningPointFrac = 0.10;   % mask bins where |dU/dbin| < this fraction of its max
SpapsTol         = 0.09;   % only used by the 'gradient' method
dQdVClamp        = 10e-8;  % only used by the 'gradient' method

% Phase-averaged results live on their own bin grid, not the sample grid
PhaseUr    = nan(Nfiles,PhaseBins);   % bin-mean Ur, V
PhaseQp    = nan(Nfiles,PhaseBins);   % bin-mean Qp, C
PhaseQ     = nan(Nfiles,PhaseBins);   % bin-mean total Q, C
PhasedQpdU = nan(Nfiles,PhaseBins);   % dQp/dU, F
PhasedQdU  = nan(Nfiles,PhaseBins);   % dQ/dU,  F
PhaseTime  = nan(Nfiles,PhaseBins);   % bin centre time within one period, s
PhaseNcyc  = zeros(1,Nfiles);         % complete cycles used per file
PhaseUg = nan(Nfiles,PhaseBins);
PhaseUd = nan(Nfiles,PhaseBins);
switch lower(ChargeDerivativeMethod)

case 'phase'
    for n = 1:Nfiles
        dt = median(diff(Time(n,:)));
        Period = 1/F_Meas_On(n);
        P = Period/dt;                     % samples per period, non-integer
        PhaseNcyc(n) = floor(Npoints/P);   % use whole cycles only

        if PhaseNcyc(n) < 1
            warning('File %d holds less than one full period; skipping dQ/dU.',n);
            continue
        end

        % Adapt the bin count so no bin can come up empty on short records
        nb = PhaseBins;
        while nb > 20 && PhaseNcyc(n)*P/nb < 4
            nb = round(nb/2);
        end
        w = min(PhaseSGWindow, 2*floor((nb-1)/4)+1);   % odd, and well under nb
        if nb ~= PhaseBins
            warning(['File %d: only %.0f samples per cycle, dropped to %d bins ' ...
                     '(SG window %d).'],n,P,nb,w);
        end

        % Full record, whole cycles only. The truncated one-period arrays are
        % deliberately NOT used here: phase averaging only pays off if it can
        % fold several cycles together.
        k  = 1:round(PhaseNcyc(n)*P);
        tt = Time(n,k);
        Ii = CurrentOn(n,k) - mean(CurrentOn(n,k));   % kill the DC offset first,
                                                      % otherwise cumtrapz ramps
        Ur = VoltageOn(n,k) - Ii*R_off(OffIndex(n));
      %   QLoss = PlasmaOffNegativeA(n) * LossTangent ...
      % * VoltageMeasured(n) ...
      % * cos(2*pi*F_Meas_On(n)*tt + Ur(n));
      QLoss = 0;
      Rp = 1 / (2 * pi() *F_Meas_On(n) *PlasmaOffNegativeA(n) * LossTangent);
      ILoss = Ur./Rp;
      QLoss = cumtrapz(tt,ILoss);
      QLoss = QLoss - mean(QLoss);
      

        Q  = cumtrapz(tt,Ii);
        Q  = Q - mean(Q);
        Qp = Q - Ur*PlasmaOffNegativeA(n) - QLoss;
        % Find Ug and Ud 
       

        % Phase bins. mod() handles the non-integer period correctly, which
        % matters here: P is 12500.19 and 7812.57 on your two files, not an
        % integer, so index-based folding would drift.
        idx = (0:numel(k)-1).';
        b   = min(floor(mod(idx/P,1)*nb)+1, nb);
        cnt = accumarray(b,1,[nb 1]);
        if any(cnt == 0)
            warning('File %d: empty phase bins, dQ/dU may be unreliable.',n);
            cnt = max(cnt,1);
        end
        Ub  = accumarray(b,Ur(:),[nb 1])./cnt;
        Qpb = accumarray(b,Qp(:),[nb 1])./cnt;
        Qb  = accumarray(b,Q(:), [nb 1])./cnt;


        % Wrap-aware Savitzky-Golay first derivative. Step size is one bin and
        % cancels in the ratio, so no scaling by dt is needed.
        [~,g] = sgolay(PhaseSGOrder,w);
        dUb  = CircDeriv(Ub, g);
        dQpb = CircDeriv(Qpb,g);

        % Mask the turning points rather than clamping them
        Cp = dQpb./dUb;
        Cp(abs(dUb) < TurningPointFrac*max(abs(dUb))) = NaN;




posMask = Ub > 0 & ~isnan(Cp);
negMask = Ub < 0 & ~isnan(Cp);

CdPosPhase(n) = median(Cp(posMask));   % robust to a few bad bins
CdNegPhase(n) = median(Cp(negMask));

% --- Invert to get Ud, then Ug, per phase bin ---
UdPhase = nan(nb,1);
UgPhase = nan(nb,1);

UdPhase(posMask) = Qpb(posMask) / CdPosPhase(n);
UdPhase(negMask) = Qpb(negMask) / CdNegPhase(n);

UgPhase(posMask) = Ub(posMask) - UdPhase(posMask);
UgPhase(negMask) = Ub(negMask) - UdPhase(negMask);

        PhaseUr(n,1:nb)    = Ub;
        PhaseQp(n,1:nb)    = Qpb;
        PhaseQ(n,1:nb)     = Qb;
        PhasedQpdU(n,1:nb) = Cp;
        PhasedQdU(n,1:nb)  = Cp + PlasmaOffNegativeA(n);   % see header note
        PhaseTime(n,1:nb)  = ((0:nb-1)+0.5)/nb*Period;
        PhaseUd(n,1:nb) = UdPhase;
        PhaseUg(n,1:nb) = UgPhase;
        % Map back onto the truncated sample grid so the Section 8 plots and the
        % Origin export keep working unchanged. Interpolation is done in phase
        % with one bin wrapped onto each end, so the seam is continuous.
        kt   = 1:TruncIndex(n);
        ph   = mod((TTrunc(n,kt) - Time(n,1))/Period, 1);
        bc   = (((0:nb-1)+0.5)/nb).';
        bcE  = [bc(end)-1; bc; bc(1)+1];
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
DischargeGraphing =0;
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
LissajousGraphing = 1;
if LissajousGraphing == 1
    for k = 40:10:Nfiles
        figure
        scatter(downsample(UrTrunc(k,1:TruncIndex(k)),20), ...
                1e6*downsample(QTrunc(k,1:TruncIndex(k)),20),20,[0 0 1]); hold on
        scatter(UrTrunc(k,point2index(k)), 1e6*QTrunc(k,point2index(k)), 'filled', 'LineWidth', 2,'Color','Green')
        scatter(UrTrunc(k,point4index(k)), 1e6*QTrunc(k,point4index(k)),'filled', 'LineWidth', 2, 'Color','Green')
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
GapVoltageGraphing = 1;
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


% Gap voltage (phase-averaged)
GapVoltageStart = min(80,Nfiles);
GapVoltageData  = [];
GapVoltageGraphing = 1;
if GapVoltageGraphing == 1
    figure(3001)
    i = 1;
    GapVoltageLegendLabels = {};
    for k = GapVoltageStart:4:Nfiles
        if all(isnan(PhaseUg(k,:)))
            continue   % file had no valid plateau Cd, e.g. no ignition detected
        end
        plot(PhaseTime(k,:), PhaseUg(k,:)); hold on
        GapVoltageLegendLabels{i} = sprintf('%.3f min',AAA(k)/60);   %#ok<SAGROW>
        GapVoltageData(2*i-1,1:PhaseBins) = PhaseTime(k,:);          %#ok<SAGROW>
        GapVoltageData(2*i,  1:PhaseBins) = PhaseUg(k,:);            %#ok<SAGROW>
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
dQdVGraphing = 1;
if dQdVGraphing == 1
    figure
    i = 1;
    dQdVLabels = {};
    for k = 1:2:Nfiles
        plot(VTrunc(k,1:TruncIndex(k))/VoltageMeasured(k), ...
             PlasmaChargeDerivitive(k,1:TruncIndex(k))*1e12); hold on
        dQdVLabels{i} = sprintf('%d V',round(VoltageMeasured(k)));   %#ok<SAGROW>
        %dQdVLabels{i} = sprintf('%d V',round(F_meas1(k)));   %#ok<SAGROW>
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
        plot(TTrunc(k,1:TruncIndex(k)),PlasmaChargeDerivitive(k,1:TruncIndex(k))*1e12); hold on
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
        PlottedData = zeros(round((Nfiles-16)/12),1000);
        for k = 16:12:Nfiles
            % plot the phase normilized to overvoltage
            
            plot(PhaseUr(k,:)/abs(max(PhaseUr(k,:))),PhasedQpdU(k,:)*1e12,'LineWidth',1.4); hold on
            phaseLabels{i} = sprintf('%d V',round(VoltageMeasured(k)));   %#ok<SAGROW>
            PlottedData(2*i-1,:) = PhaseUr(k,:)/abs(max(PhaseUr(k,:)));
            PlottedData(2*i,:) = PhasedQpdU(k,:)*1e12;
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
DielectricCapacitancePosTrip  = zeros(1,NTrip);
DielectricCapacitanceTripPosSTD = zeros(1,NTrip);
DielectricCapacitanceNegTrip  = zeros(1,NTrip);
DielectricCapacitanceNegTripSTD = zeros(1,NTrip);
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
    DielectricCapacitanceTripPos(n)  = mean(PlasmaOnPositiveA(idx));
    DielectricCapacitanceTripPosSTD(n)= std(PlasmaOnPositiveA(idx));
    DielectricCapacitanceTripNeg(n)  = mean(PlasmaOnNegativeA(idx));
    DielectricCapacitanceTripNegSTD(n)= std(PlasmaOnNegativeA(idx));
    % FIX: guarded. The original referenced PlasmaExtentPositive unconditionally
    % while the block that defines it was commented out, so this loop always
    % errored with "Unrecognized function or variable 'PlasmaExtentPositive'".
    if ComputePlasmaExtent == 1
        PlasmaExtentPositiveTrip(n)   = mean(PlasmaExtentPositive(idx));
        PlasmaExtentPositiveTripSTD(n)= std(PlasmaExtentPositive(idx));
    end
end
OriginStuff = [VoltageMeasuredTrip;VoltageMeasuredTripSTD;DielectricCapacitanceTripPos;DielectricCapacitanceTripPosSTD;DielectricCapacitanceTripNeg;DielectricCapacitanceTripNegSTD]';
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

figure;
scatter(UrMeasured, (UrMeasured + IgnitionVoltagePositive)/2  ); hold on
scatter(UrMeasured,(UrMeasured - IgnitionVoltageNegative)/2 ); hold on
legend('positive','negative')
title('Ignition Voltage Difference')
figure;
scatter(UrMeasured, IgnitionVoltagePositive ); hold on
scatter(UrMeasured,IgnitionVoltageNegative ); hold on
legend('positive','negative')
title('Ignition Voltages')


% figure;
% scatter(UrMeasured, 0.5*(1-alpha).*(UrMeasured + IgnitionVoltagePositive) ); hold on
% scatter(UrMeasured,0.5*(1-alpha).* (UrMeasured - IgnitionVoltageNegative) ); hold on
% legend('positive','negative')
% title('Callibration')


%% Section 11 - Export for Origin
% FIX: the original ran this unconditionally using an index vector called Index
% that is never created anywhere in the script, so the last six lines always
% threw "Unrecognized function or variable 'Index'". Now gated, and Index is
% built explicitly. Set the file numbers you actually want to export.
ExportForOrigin = 1;
Index = [61,67,75,83,88,91];   % e.g. Index = [1 20 40 60 80 100];

if ExportForOrigin == 1
    if numel(Index) < 4
        error('Index must contain at least four file numbers to export.');
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

function [fs,V_mag,I_mag] = VoltageTransform(t,V,I)
t = time;        % Time vector (s)
V = voltage;     % Voltage signal
I = current;     % Current signal

%% Sampling parameters
Fs = 1/mean(diff(t));      % Sampling frequency (Hz)
N  = length(t);

%% Remove DC offset
V = V - mean(V);
I = I - mean(I);

%% FFT
V_fft = fft(V);
I_fft = fft(I);

%% Single-sided spectrum
f = Fs*(0:N/2)/N;

V_mag = abs(V_fft/N);
I_mag = abs(I_fft/N);

V_mag = V_mag(1:N/2+1);
I_mag = I_mag(1:N/2+1);

% Double non-DC components
V_mag(2:end-1) = 2*V_mag(2:end-1);
I_mag(2:end-1) = 2*I_mag(2:end-1);

%% Plot
figure;

subplot(2,1,1)
plot(f,V_mag,'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel('|V(f)|')
title('Voltage FFT')
grid on
xlim([0 Fs/2])

subplot(2,1,2)
plot(f,I_mag,'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel('|I(f)|')
title('Current FFT')
grid on
xlim([0 Fs/2])
end



Urint = interp1(T3,URR,T1);

