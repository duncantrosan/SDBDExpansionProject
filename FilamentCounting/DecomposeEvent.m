function [Comp,Fit] = DecomposeEvent(y,Tp,Sigma,Opt)
% Split one discharge event into individual current pulses.
%   y:     the event's discharge current (row, A), polarity flipped so pulses are positive,
%          local baseline removed; sample k is at time (k-1)*Tp.dt.
%   Tp:    unit-charge pulse template from BuildPulseTemplate.
%   Sigma: noise sigma of the current (A).
% Model: y(t) = sum_j q_j h(t - tau_j) + c with q_j >= 0. Pulses are added one at a time
% where the matched-filter gain (drop in squared error) is largest. A candidate is accepted
% only if
%   - the gain exceeds (KDet*Sigma)^2, i.e. a KDet-sigma matched-filter detection,
%   - its peak current exceeds Eta times the already-fitted waveform around it (+-10 ns), and
%   - within ShadowTime after an accepted pulse, its peak exceeds Shadow times that pulse's.
% The last two tests stop the fit from explaining the template's shape error with extra
% pulses: large events ring longer and at a different period than the small pulses the
% template is built from, by up to a third of their peak, and without these tests every
% ringing lobe would become a "pulse". A real pulse hidden under that ringing is lost too.
% After each addition all charges are re-fitted (non-negative least squares) and the new
% pulse time is refined. Pulses closer than MinSep are not allowed.
% Returns Comp (table: tau_ns from segment start, q_nC) and Fit (fitted waveform).
Def = struct('KDet',6,'Eta',0.5,'Shadow',0.35,'ShadowTime',400e-9,'MinSep',4e-9,'MaxComp',40);
if nargin < 4, Opt = struct; end
f = fieldnames(Def);
for k = 1:numel(f), if ~isfield(Opt,f{k}), Opt.(f{k}) = Def.(f{k}); end, end
dt = Tp.dt; Up = Tp.Up;
y = y(:); L = numel(y); tk = (0:L-1)'*dt;
Pre = ceil(-Tp.tFine(1)/dt); Post = floor(Tp.tFine(end)/dt);
% coarse templates for each sub-sample phase: G{p}(i) = h(i*dt - p*dt/Up), i = -Pre..Post
G = cell(1,Up);
for p = 0:Up-1, G{p+1} = interp1(Tp.tFine,Tp.h,((-Pre:Post) - p/Up)*dt,'linear',0)'; end

Tau = zeros(0,1); q = zeros(0,1);
Res = y - mean(y);
HPk = max(Tp.h);
Env = zeros(L,1);                                       % |fit| envelope, +-10 ns
Sh = zeros(L,1);                                        % shadow of earlier pulses
for It = 1:Opt.MaxComp
    Best = -Inf; BestTau = NaN;
    for p = 0:Up-1
        g = G{p+1};
        c = conv(Res, flipud(g));                       % c(m+Post+1) = sum_k Res(k) g(k-1-m)
        e = conv(ones(L,1), flipud(g.^2));
        m = (0:L-1)';                                   % peak at sample m (time m*dt + p*dt/Up)
        cm = c(m + Post + 1); em = e(m + Post + 1);
        Gain = (max(cm,0).^2)./max(em,eps);
        Gain(cm./max(em,eps)*HPk < max(Opt.Eta*Env,Sh)) = 0;   % would only patch the template's shape error
        tc = (m + p/Up)*dt;
        if ~isempty(Tau)
            Gain(min(abs(tc - Tau'),[],2) < Opt.MinSep) = 0;
        end
        [gm,im] = max(Gain);
        if gm > Best, Best = gm; BestTau = tc(im); end
    end
    if Best < (Opt.KDet*Sigma)^2, break, end
    Tau(end+1,1) = BestTau; %#ok<AGROW>
    [q,Res,Tau] = Refit(y,tk,Tp,Tau,Opt.MinSep,numel(Tau));
    Env = movmax(abs(y - Res - mean(y - Res)), 2*round(10e-9/dt)+1);
    Sh = zeros(L,1);
    for j = 1:numel(Tau)
        In = tk > Tau(j) & tk <= Tau(j) + Opt.ShadowTime;
        Sh(In) = max(Sh(In), Opt.Shadow*q(j)*HPk);
    end
end
if ~isempty(Tau), [q,Res,Tau] = Refit(y,tk,Tp,Tau,Opt.MinSep,1:numel(Tau)); end
[Tau,o] = sort(Tau); q = q(o);
Comp = table(Tau*1e9, q*1e9, 'VariableNames', {'tau_ns','q_nC'});
Fit = y - Res;
end

function [q,Res,Tau] = Refit(y,tk,Tp,Tau,MinSep,Which)
% Non-negative charges plus a free constant; then walk the times of the pulses in Which
% one fine step at a time while the squared error keeps dropping.
Step = Tp.dt/Tp.Up;
[q,Rss] = Solve(y,tk,Tp,Tau);
for j = Which(:)'
    for d = [-1 1]
        while true
            T2 = Tau; T2(j) = T2(j) + d*Step;
            Others = T2([1:j-1 j+1:end]);
            if any(abs(Others - T2(j)) < MinSep), break, end
            [q2,R2] = Solve(y,tk,Tp,T2);
            if R2 >= Rss, break, end
            Tau = T2; q = q2; Rss = R2;
        end
    end
end
keep = q > 0;                                           % drop pulses NNLS switched off
Tau = Tau(keep);
[q,~,Res] = Solve(y,tk,Tp,Tau);
end

function [q,Rss,Res] = Solve(y,tk,Tp,Tau)
H = interp1(Tp.tFine,Tp.h,tk - Tau','linear',0);
A = [H ones(numel(y),1) -ones(numel(y),1)];
Scale = sqrt(sum(A.^2,1)); Scale(Scale == 0) = 1;      % unit columns keep NNLS well conditioned
x = lsqnonneg(A./Scale,y,optimset('Display','off'))./Scale';
q = x(1:numel(Tau));
Res = y - A*x;
Rss = Res'*Res;
end
