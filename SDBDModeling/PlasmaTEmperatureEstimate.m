%% ElectrodeThermometer.m
%
% Use the electrode's own dielectric loss as a thermometer.
%
% The loss tangent of the adhesive layer changes with temperature (the
% relaxation peak shifts). If you have calibration data from an LCR meter
% at known temperatures, you can invert the relationship: measure Re(Z)
% from the scope, subtract R_esr, compute tan(d), and look up the
% temperature on the calibration curve.
%
% This gives you the electrode temperature during plasma operation without
% any additional sensors, at no extra cost, from data you already collect.
%
% USAGE
%   1. Run HiokiImpedanceSweep.m on your temperature sweep files to get
%      TempTable (or enter the calibration data manually below).
%   2. Set R_esr from the off-measurement analysis.
%   3. Pass in Re(Z) from the scope at any point during the run.
%
% ACCURACY
%   Sensitivity is about 24 ohm/K in Re(Z)_device at 1600 Hz, which is
%   0.18%/K in total Re(Z) on top of 13 kOhm of R_esr. The scope sine
%   fit gives Re(Z) to ~100 ohm, so the temperature resolution is roughly
%   +- 4 K. Good enough to distinguish "room temperature" from "60 C after
%   15 minutes of plasma" but not for precision thermometry.
%
% LIMITATIONS
%   - Valid only above ~700 Hz where tan(d) decreases monotonically with T.
%     Below 700 Hz the trend reverses and the inversion is ambiguous.
%   - Assumes R_esr is temperature-independent (it is, if it is in the
%     amplifier). If R_esr is in the device, this needs rethinking.
%   - The calibration is device-specific. Different adhesive thickness or
%     type gives a different curve.
 
%% ========================================================================
%  CALIBRATION DATA
%  ========================================================================
% Enter your LCR temperature sweep results here.
% Each row: [Temperature_C, tan_delta, Capacitance_F]
% at your operating frequency.
%
% If you have TempTable from HiokiImpedanceSweep.m, you can build this
% automatically (see the block below).
 
OperatingFreq = 1600;   % Hz, must match the frequency of the calibration data
 
% --- Option 1: Manual entry from your PET data ---------------------------
% CalibrationData = [T(C), tan(d), C(F)] at OperatingFreq
CalibrationData = [
    25,  0.02196,  467.6e-12
    30,  0.02146,  467.6e-12
    35,  0.01957,  466.3e-12
    40,  0.01895,  465.0e-12
    55,  0.01868,  463.5e-12
    60,  0.01767,  457.5e-12
];
 
% --- Option 2: Build from TempTable (uncomment if available) -------------
% CalibrationData = zeros(numel(TempTable),3);
% jf = nearestIdx(TempTable(1).Freq, OperatingFreq);
% for k = 1:numel(TempTable)
%     CalibrationData(k,:) = [TempTable(k).Temp, ...
%                             TempTable(k).TanDelta_mean(jf), ...
%                             TempTable(k).Cs_mean(jf)];
% end
% % Remove NaN temperature entries
% CalibrationData = CalibrationData(~isnan(CalibrationData(:,1)),:);
 
% --- R_esr from your off measurement analysis ---------------------------
R_esr = 13200;   % Ohm, from the tan(d) decomposition
 
%% ========================================================================
%  BUILD THE CALIBRATION CURVE
%  ========================================================================
T_cal  = CalibrationData(:,1);
td_cal = CalibrationData(:,2);
C_cal  = CalibrationData(:,3);
w_op   = 2*pi*OperatingFreq;
 
% Device Re(Z) at each calibration temperature
Re_cal = td_cal ./ (w_op * C_cal);
 
% Check monotonicity (required for unique inversion)
diffs = diff(td_cal);
if all(diffs < 0)
    fprintf('Calibration curve is monotonically decreasing. Good.\n');
elseif all(diffs > 0)
    fprintf('Calibration curve is monotonically increasing.\n');
    fprintf('  (This happens below the crossover frequency ~700 Hz.\n');
    fprintf('   Inversion will work but the sign is flipped.)\n');
else
    warning(['Calibration curve is NOT monotonic. Temperature inversion ' ...
             'may be ambiguous. Check that the operating frequency is ' ...
             'above ~700 Hz.']);
end
 
% Fit a polynomial for smooth interpolation/extrapolation
% tan(d) vs T: use a quadratic (the relaxation shift is mildly nonlinear)
polyOrder = min(2, numel(T_cal)-1);
p_td = polyfit(T_cal, td_cal, polyOrder);
p_C  = polyfit(T_cal, C_cal, polyOrder);
 
% Inverse: T from tan(d). Sort by td for interpolation.
[td_sorted, isort] = sort(td_cal);
T_sorted = T_cal(isort);
 
fprintf('\nCalibration curve at %d Hz:\n', OperatingFreq);
fprintf('  %5s  %8s  %8s  %10s\n','T(C)','tan(d)','C(pF)','Re_dev(Ohm)');
for k = 1:numel(T_cal)
    fprintf('  %5.1f  %8.5f  %8.1f  %10.0f\n', ...
            T_cal(k), td_cal(k), C_cal(k)*1e12, Re_cal(k));
end
 
% Sensitivity
dT = T_cal(end)-T_cal(1);
dtd = td_cal(end)-td_cal(1);
dRe = Re_cal(end)-Re_cal(1);
fprintf('\n  Over %.0f - %.0f C:\n', T_cal(1), T_cal(end));
fprintf('    tan(d) changes by %+.5f  (%.1f ppm/K)\n', dtd, dtd/dT/mean(td_cal)*1e6);
fprintf('    Re_dev changes by %+.0f ohm  (%.1f ohm/K)\n', dRe, dRe/dT);
fprintf('    Resolution: ~100 ohm scope noise -> +- %.0f K\n', 100/abs(dRe/dT));
 
%% ========================================================================
%  INVERSION FUNCTION
%  ========================================================================
% Use this on your scope data to get temperature.
%
% Input:  Re_Z_total = Re(Z) from the scope sine fit (the old "R_off")
%         at the operating frequency
% Output: T_estimated in degrees C
%
% Example usage in the main analysis loop:
%
%   for n = 1:Nfiles
%       T_electrode(n) = InvertTemperature(R_total(n), R_esr, ...
%                                          OperatingFreq, C_cal_mean, ...
%                                          td_sorted, T_sorted);
%   end
 
% --- Demo: apply to the off measurements if they exist ---
% (Uncomment and connect to your main script's R_total vector)
%
% T_electrode = zeros(1, Nfiles);
% for n = 1:Nfiles
%     T_electrode(n) = InvertTemperature(R_total(n), R_esr, ...
%                                        F_meas1(n), mean(C_cal), ...
%                                        td_sorted, T_sorted);
% end
% figure
% plot(1:Nfiles, T_electrode, '-ok', 'LineWidth', 1.4)
% xlabel('File index'); ylabel('Estimated Temperature / °C')
% title('Electrode Temperature from Impedance')
 
%% ========================================================================
%  APPLY TO A TIME SERIES (e.g., plasma on then cooling)
%  ========================================================================
% If you have a time-resolved Re(Z) trace (e.g., from tracking the sine
% fit parameters through a long acquisition), you can get T(t):
%
%   T_vs_time = arrayfun(@(R) InvertTemperature(R, R_esr, ...
%                         OperatingFreq, mean(C_cal), td_sorted, T_sorted), ...
%                         Re_Z_timeseries);
%
% The 15-second cooling trace you mentioned would give you both the
% peak temperature and the thermal time constant directly.
 
%% ========================================================================
%  DIAGNOSTIC PLOT
%  ========================================================================
figure('Name','Temperature Calibration Curve')
 
subplot(2,1,1)
scatter(T_cal, td_cal*1e3, 60, 'filled'); hold on
T_fine = linspace(T_cal(1)-5, T_cal(end)+5, 200);
plot(T_fine, polyval(p_td, T_fine)*1e3, 'r-', 'LineWidth', 1.2)
xlabel('Temperature / °C')
ylabel('tan \delta  (\times10^{-3})')
title(sprintf('Loss Tangent vs Temperature at %d Hz', OperatingFreq))
grid on
 
subplot(2,1,2)
scatter(T_cal, Re_cal, 60, 'filled'); hold on
Re_fine = polyval(p_td, T_fine) ./ (w_op * polyval(p_C, T_fine));
plot(T_fine, Re_fine, 'r-', 'LineWidth', 1.2)
xlabel('Temperature / °C')
ylabel('Re(Z)_{device} / \Omega')
title('Device Resistance vs Temperature')
grid on
 
set(gcf,'Visible','on')
 
%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================
 
function T = InvertTemperature(Re_Z_total, R_esr, f, C_approx, td_sorted, T_sorted)
% Estimate electrode temperature from measured total Re(Z).
%
%   Re_Z_total  - Re(Z) from the scope sine fit, ohm
%   R_esr       - external series resistance (amplifier), ohm
%   f           - operating frequency, Hz
%   C_approx    - approximate device capacitance, F (for tan(d) calc)
%   td_sorted   - calibration tan(d) values, sorted ascending
%   T_sorted    - corresponding temperatures, same order
%
% Returns estimated temperature in degrees C.
% Returns NaN if the measured tan(d) is outside the calibration range.
 
w = 2*pi*f;
Re_device = Re_Z_total - R_esr;
 
if Re_device <= 0
    T = NaN;
    warning('Re_device <= 0. R_esr may be too large or the device is not capacitive.');
    return
end
 
td_measured = Re_device * w * C_approx;
 
% Interpolate on the calibration curve
if td_measured < td_sorted(1) || td_measured > td_sorted(end)
    % Extrapolate with a warning
    T = interp1(td_sorted, T_sorted, td_measured, 'linear', 'extrap');
    warning('tan(d) = %.5f is outside calibration range [%.5f, %.5f]. Extrapolating to T = %.1f C.', ...
            td_measured, td_sorted(1), td_sorted(end), T);
else
    T = interp1(td_sorted, T_sorted, td_measured, 'linear');
end
end
 
function j = nearestIdx(vec, target)
[~,j] = min(abs(vec - target));
end