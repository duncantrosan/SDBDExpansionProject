function Fq = FitChargeQuanta(q,QRange,NMax,NBoot)
% Test whether pulse charges come in multiples of a single-filament charge q1.
%   q:      pulse charges (nC, positive).
%   QRange: [qCut qMax] nC. Only charges in this range are fitted (both models are truncated
%           there): detection is incomplete below qCut, and above qMax the multiplicities
%           are too high for humps to be resolvable anyway.
%   NMax:   largest multiplicity considered (default 8).
%   NBoot:  parametric-bootstrap replicates for the p-value (default 100; 0 skips it).
% Models
%   Smooth - two-component lognormal mixture (5 parameters): flexible, but no humps.
%   Quanta - sum_n w_n N(n*q1, s0^2 + n*s1^2), n = 1..NMax, zero-truncated Poisson weights
%            (4 parameters: q1, s0, s1, lambda). s0 is the measurement spread, s1 the
%            filament-to-filament spread, which adds up over n filaments. Both are held
%            below q1/2 so that the humps stay resolvable: otherwise the "quantized" model
%            can turn itself into a smooth one and the test means nothing.
% Returns Fq with the fitted parameters, log-likelihoods, BIC of both models, the bootstrap
% p-value of the likelihood gain under the Smooth model, and for each input charge the
% most probable multiplicity N and its posterior probability PN (NaN outside QRange).
if nargin < 3 || isempty(NMax), NMax = 8; end
if nargin < 4 || isempty(NBoot), NBoot = 100; end
q = q(:); In = q >= QRange(1) & q <= QRange(2); x = q(In);
[Ps,Ls] = FitSmooth(x,QRange);
[Pq,Lq] = FitQuanta(x,QRange,NMax);
n = numel(x);
Fq.q1 = Pq(1); Fq.s0 = Pq(2); Fq.s1 = Pq(3); Fq.Lambda = Pq(4);
Fq.LogLSmooth = Ls; Fq.LogLQuanta = Lq;
Fq.BICSmooth = -2*Ls + 5*log(n); Fq.BICQuanta = -2*Lq + 4*log(n);
Fq.NFit = n; Fq.QRange = QRange; Fq.NMax = NMax;
Fq.ParSmooth = Ps; Fq.ParQuanta = Pq;
Fq.PBoot = NaN; Fq.GainBoot = [];
if NBoot > 0
    Gain = zeros(NBoot,1);
    for b = 1:NBoot
        xb = SampleSmooth(Ps,QRange,n);
        [~,Lsb] = FitSmooth(xb,QRange);
        [~,Lqb] = FitQuanta(xb,QRange,NMax);
        Gain(b) = Lqb - Lsb;
    end
    Fq.GainBoot = Gain;
    Fq.PBoot = (1 + sum(Gain >= Lq - Ls))/(NBoot + 1);
end
[~,Comp] = QuantaPdf(q,Pq,NMax);
[PN,N] = max(Comp./sum(Comp,2),[],2);
N(~In) = NaN; PN(~In) = NaN;
Fq.N = N; Fq.PN = PN;
end

function [P,L] = FitSmooth(x,R)
lx = log(x); m = mean(lx); s = std(lx);
Starts = [m-0.5*s 0.5*s m+0.7*s 0.6*s 0.5; m-0.2*s 0.7*s m+1.2*s 0.4*s 0.7; m 0.8*s m 1.5*s 0.5];
Best = Inf;
for k = 1:size(Starts,1)
    z0 = [Starts(k,1) log(Starts(k,2)) Starts(k,3) log(Starts(k,4)) log(Starts(k,5)/(1-Starts(k,5)))];
    [z,f] = fminsearch(@(z) -SmoothLL(x,R,z),z0,optimset('MaxFunEvals',4000,'MaxIter',4000,'Display','off'));
    if f < Best, Best = f; Zb = z; end
end
P = Zb; L = -Best;
end

function L = SmoothLL(x,R,z)
w = 1/(1+exp(-z(5))); s1 = exp(z(2)); s2 = exp(z(4));
f = w*lognpdf(x,z(1),s1) + (1-w)*lognpdf(x,z(3),s2);
In = w*(logncdf(R(2),z(1),s1) - logncdf(R(1),z(1),s1)) + (1-w)*(logncdf(R(2),z(3),s2) - logncdf(R(1),z(3),s2));
L = sum(log(max(f,1e-300))) - numel(x)*log(max(In,1e-300));
end

function x = SampleSmooth(z,R,n)
% Exact draws from the mixture truncated to R (inverse CDF per component, no rejection,
% which would never finish if the fit puts most of its mass outside R).
Mu = [z(1) z(3)]; Sd = exp([z(2) z(4)]);
Lo = logncdf(R(1),Mu,Sd); Hi = logncdf(R(2),Mu,Sd);
w = 1/(1+exp(-z(5))); W = [w 1-w].*(Hi - Lo);
c = 1 + (rand(n,1) > W(1)/sum(W));
u = Lo(c)' + rand(n,1).*(Hi(c) - Lo(c))';
x = logninv(u,Mu(c)',Sd(c)');
x = min(max(x,R(1)),R(2));
end

function [P,L] = FitQuanta(x,R,NMax)
Best = Inf;
for q1 = linspace(0.3,2.5,12)
    z0 = [log(q1) log(0.1*q1) log(0.15*q1) log(1.5)];
    [z,f] = fminsearch(@(z) -QuantaLL(x,R,z,NMax),z0,optimset('MaxFunEvals',3000,'MaxIter',3000,'Display','off'));
    if f < Best, Best = f; Zb = z; end
end
P = exp(Zb); L = -Best;
end

function L = QuantaLL(x,R,z,NMax)
P = exp(z);
if P(1) < 0.3 || P(1) > 4 || P(2) < 0.01 || P(2) > P(1)/2 || P(3) > P(1)/2 || P(4) > 50
    L = -1e10; return
end
[f,~,In] = QuantaPdf(x,P,NMax,R);
L = sum(log(max(f,1e-300))) - numel(x)*log(max(In,1e-300));
end

function [f,Comp,In] = QuantaPdf(x,P,NMax,R)
nn = 1:NMax;
w = P(4).^nn./factorial(nn); w = w/sum(w);
mu = nn*P(1); sd = sqrt(P(2)^2 + nn*P(3)^2);
Comp = w.*normpdf(x(:),mu,sd);
f = sum(Comp,2);
if nargin > 3, In = sum(w.*(normcdf(R(2),mu,sd) - normcdf(R(1),mu,sd))); else, In = 1; end
end
