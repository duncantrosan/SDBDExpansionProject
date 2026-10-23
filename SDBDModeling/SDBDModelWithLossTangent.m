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
ThreshInput       = 0.00035;   % A (0.35 mA)
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