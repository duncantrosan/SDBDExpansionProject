Res = real(A(2,:));
Reac = -imag(A(2,:));
Freq = A(1,:);




% Convert frequency to angular frequency
w = 2 * pi * Freq;

% Define residual function for fitting both real and imaginary parts
residualFunction = @(params) [ ...
    params(2) * w.^2 .* (params(1)*1e-12).^2 ./ (w.^2 .* (params(1)*1e-12).^2 + params(2).^2) - Res; ...
    w .* (params(1)*1e-12) .* params(2).^2 ./ (w.^2 .* (params(1)*1e-12).^2 + params(2).^2) - 1 ./ (w .* (params(3)*1e-12)) - Reac ...
];

% Weighting factor to balance real/imag part residuals
wgt_real = 1;
wgt_imag = 1;

residualFunction = @(params) [ ...
    wgt_real * (params(2) * w.^2 .* params(1).^2 ./ (w.^2 .* params(1).^2 + params(2).^2) - Res); ...
    wgt_imag * (w .* params(1) .* params(2).^2 ./ (w.^2 .* params(1).^2 + params(2).^2) - 1 ./ (w .* params(3)) - Reac) ...
];

% Initial guesses in pF, R in Ohm
initialGuess = [2000, 25000, 600];  % means C1=2000pF=2nF, R=5000 Ohm, C2=1500pF=1.5nF
lb = [0.0001, 0.011, 0.01];       % 1 pF minimum
ub = [1e5, 1e10, 1e5]; % 100 nF maximum


% Options for solver display and tolerance
options = optimoptions('lsqnonlin', 'Display', 'iter', 'TolFun',1e-12, 'TolX',1e-12);

% Perform nonlinear least squares fitting
fittedParams = lsqnonlin(residualFunction, initialGuess, lb, ub, options);

% Extract fitted parameters
C1_fit = fittedParams(1);
R_fit = fittedParams(2);
C2_fit = fittedParams(3);

% Display fitted results
fprintf('Fitted C1 = %.3f pF\n', C1_fit);
fprintf('Fitted R = %.3f kOhm\n', R_fit/1000);
fprintf('Fitted C2 = %.3f pF\n', C2_fit);


% func = @(x) a * x + b;
% 
% func1 = @(w) w^2*C1^2*R/(w^2 * C1^2 + R^2);
% func2 = @(w) w*C1 *R^2 / (w^2 * C1^2 + R^2) - 1/(w*C2); 
% 
% w = 2*pi();
% Res = func1;
% Reac = func2;


%% Plotting 
% Convert fitted capacitances back to Farads
C1 = fittedParams(1) * 1e-12;
R = fittedParams(2);
C2 = fittedParams(3) * 1e-12;

% Compute modeled real and imaginary parts using the fitted parameters
Res_fit = R .* w.^2 .* C1.^2 ./ (w.^2 .* C1.^2 + R.^2);
Reac_fit = w .* C1 .* R.^2 ./ (w.^2 .* C1.^2 + R.^2) - 1 ./ (w .* C2);

figure;

% Real part plot
subplot(2,1,1);
plot(Freq, Res, 'bo', 'DisplayName', 'Experimental Real');
hold on;
plot(Freq, Res_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Real');
xlabel('Frequency (Hz)');
ylabel('Real Part of Impedance (Ohm)');
legend;
grid on;
title('Real Part: Experimental vs. Fitted');

% Imaginary part plot
subplot(2,1,2);
plot(Freq, Reac, 'bo', 'DisplayName', 'Experimental Imag');
hold on;
plot(Freq, Reac_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Imag');
xlabel('Frequency (Hz)');
ylabel('Imaginary Part of Impedance (Ohm)');
legend;
grid on;
title('Imaginary Part: Experimental vs. Fitted');


% Define circuit parameters (adjust these to experiment with fitting)
C1 = 10e-12;  % 2 nF
R = 25000;   % 5 kohm
C2 = 600e-12; % 1.5 nF

% Angular frequency
w = 2 * pi * linspace(0,100,100);

% Theoretical impedance from circuit model (real and imaginary parts)
Res_fit = R .* w.^2 .* C1.^2 ./ (w.^2 .* C1.^2 + R.^2);
Reac_fit = w .* C1 .* R.^2 ./ (w.^2 .* C1.^2 + R.^2) - 1 ./ (w .* C2);

% Plot experimental data and fitted curve
figure;

% Real Part
subplot(2,1,1);
plot(w/(2*pi()), Res_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Real');
xlabel('Frequency (Hz)');
ylabel('Real Part of Impedance (Ohm)');
title('Real Part: Experimental vs. Fitted');
legend show;
grid on;

% Imaginary Part
subplot(2,1,2);
plot(w/(2*pi()), Reac_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Imaginary');
xlabel('Frequency (Hz)');
ylabel('Imaginary Part of Impedance (Ohm)');
title('Imaginary Part: Experimental vs. Fitted');
legend show;
grid on;