function Tp = BuildPulseTemplate(Segs,dt,Up)
% Average single-pulse response of the current measurement, normalised to unit charge.
%   Segs: one pulse per row (polarity flipped so the pulse is positive, local baseline
%         removed), roughly aligned on their peaks.
%   Up:   oversampling factor for sub-sample alignment (default 8).
% Each pulse is spline-upsampled, scaled to the current mean by least squares (robust for
% small, noisy pulses), and re-aligned to the mean by cross-correlation; three passes.
% Returns Tp.h (A per C, fine grid), Tp.tFine (s, 0 at the peak), Tp.dt, Tp.Up, Tp.N,
% Tp.Spread (rms deviation of the aligned, scaled pulses from the mean, relative to its peak).
if nargin < 3, Up = 8; end
[n,L] = size(Segs);
x = 0:L-1; xf = 0:1/Up:L-1;
F = zeros(n,numel(xf));
for k = 1:n, F(k,:) = interp1(x,Segs(k,:),xf,'spline'); end
Ref = mean(F,1);
MaxLag = 4*Up;
for Pass = 1:3
    for k = 1:n
        [c,lag] = xcorr(F(k,:),Ref,MaxLag);
        [~,j] = max(c);
        F(k,:) = ShiftZero(F(k,:),-lag(j));
        F(k,:) = F(k,:) * (F(k,:)*Ref')/(F(k,:)*F(k,:)') ;   % scale pulse onto the mean
    end
    Ref = mean(F,1);
end
[~,ip] = max(Ref);
Tp.tFine = ((0:numel(xf)-1) - (ip-1))*dt/Up;
Tp.h = Ref / (trapz(Ref)*dt/Up);
Tp.dt = dt; Tp.Up = Up; Tp.N = n;
Tp.Spread = sqrt(mean((F - Ref).^2,'all'))/max(Ref);
end

function y = ShiftZero(x,s)
% Shift a row by s samples (positive = later), padding with zeros.
y = zeros(size(x));
if s >= 0, y(1+s:end) = x(1:end-s); else, y(1:end+s) = x(1-s:end); end
end
