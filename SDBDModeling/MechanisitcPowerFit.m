function [a, V0, R2] = MechanisticPowerFit(V, P, Cd_slope, Cd_intercept)
% Fit measured DBD power to the mechanistic form  P = a*(V - V0)*V^2
% derived from linear Cd(V) and capacitive energy scaling E ~ Cd(V)*V^2.
%
%   V, P          : column vectors of drive amplitude [V] and cycle power [W]
%   Cd_slope      : (optional) linear Cd slope k  [F/V], e.g. 0.235e-12
%   Cd_intercept  : (optional) V0 from the Cd(V) plot [V], e.g. 1172
%                   If given, printed for comparison with the fitted V0.
%
% Returns fitted prefactor a, threshold V0, and R^2.
 
V = V(:); P = P(:);
 
% Two-parameter mechanistic model
model = @(p, V) p(1) .* (V - p(2)) .* V.^2;      % p = [a, V0]
 
% Starting guesses: V0 just below the lowest measured voltage; a from a data point
p0  = [P(end)/((V(end)-0.9*min(V))*V(end)^2), 0.9*min(V)];
opt = optimset('Display','off');
p   = lsqcurvefit(model, p0, V, P, [], [], opt);
 
a  = p(1);
V0 = p(2);
res = P - model(p, V);
R2  = 1 - sum(res.^2)/sum((P - mean(P)).^2);
 
fprintf('P = a*(V - V0)*V^2   ->   a = %.4g,  V0 = %.1f V,  R^2 = %.5f\n', a, V0, R2);
if nargin >= 4
    fprintf('Cd(V) intercept from Fig.7: V0 = %.1f V   (agreement check)\n', Cd_intercept);
end
 
% Plot data + fit
Vfit = linspace(min(V), max(V), 200)';
figure; plot(V, P, 'o', Vfit, model(p, Vfit), '-', 'LineWidth', 1.2);
xlabel('Drive amplitude V'); ylabel('Cycle power P');
legend('Measured','a(V-V_0)V^2','Location','northwest'); grid on;
end