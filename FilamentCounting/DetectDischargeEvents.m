function Ev = DetectDischargeEvents(t,V,I,Ccell,Opt)
% Find discharge current pulses in one capture.
%   t (s), V (V, already scaled by the HV probe ratio), I (A), Ccell (F) from the cell
%   capacitance check. Opt overrides any of the defaults below.
%
% Steps
%   1. Displacement current Ccell*dV/dt and the CH2 DC offset are removed (offset = median
%      of the samples before the supply switches and of the last 20% of the record).
%   2. Noise sigma is a robust (MAD) estimate on those same quiet samples; the detection
%      threshold is KSigma*sigma, so it follows the V/div setting instead of a fixed 20 mA.
%   3. Masked: the supply switch-on burst (zero-mean ~60 MHz ringing while V is still ~0) up
%      to the moment |V| first exceeds EMIEndVoltage, and EdgeGuard at both record ends.
%   4. Samples above threshold are grouped into events separated by quiet gaps > GapTime. A
%      group starting less than MergeTime after the previous one, with a peak below
%      MergeFrac of that one's, is a late ringing lobe of it and is merged into it.
%      Each event is integrated from PreTime before its first to PostTime after its last
%      above-threshold sample, against a local baseline taken over BaseTime before it, so
%      the slow (~1 mA) baseline wander of CH2 does not leak into the charge.
%   5. Expected polarity is the sign of dV/dt (positive pulses while V rises).
%
% Returns struct Ev with Ic (discharge current, A), Offset, Sigma, Threshold (A), tSwitch (s),
% Masked (logical), and table Ev.Events: Start, Stop, Peak (sample indices), tPeak_us (after
% tSwitch), V_kV at the event start, Pol (+1/-1), Ipk_A (signed), Q_nC (signed, apparent
% charge), Base_A and Width_ns (FWHM of the largest pulse).
Def = struct('KSigma',6,'GapTime',100e-9,'MergeTime',300e-9,'MergeFrac',0.35,'PreTime',10e-9,'PostTime',100e-9, ...
    'BaseTime',200e-9,'EMIBefore',2e-6,'EMIEndVoltage',2500,'EdgeGuard',20e-6,'SmoothTime',1.6e-6);
if nargin < 5, Opt = struct; end
f = fieldnames(Def);
for k = 1:numel(f), if ~isfield(Opt,f{k}), Opt.(f{k}) = Def.(f{k}); end, end

dt = t(2)-t(1);
ns = @(x) max(1,round(x/dt));
N = numel(t);
dVdt = smooth(gradient(V)./dt, 2*floor(ns(Opt.SmoothTime)/2)+1, 'sgolay', 3)';

% Supply switch: first departure of V from its level at the start of the record
Vpre = median(V(t < t(1) + 20e-6));
iSw = find(abs(V - Vpre) > 300, 1);
if isempty(iSw), iSw = round(N/2); end
iHigh = iSw - 1 + find(abs(V(iSw:end) - Vpre) > Opt.EMIEndVoltage, 1);
if isempty(iHigh), iHigh = N; end                   % never high enough to ignite: mask it all
Masked = false(1,N);
Masked(max(1,iSw-ns(Opt.EMIBefore)):max(iHigh,iSw+ns(3e-6))) = true;
Masked([1:ns(Opt.EdgeGuard) N-ns(Opt.EdgeGuard)+1:N]) = true;

Quiet = false(1,N);
Quiet(1:max(1,iSw-ns(5e-6))) = true;
Quiet(round(0.8*N):N) = true;
Offset = median(I(Quiet));
Ic = I - Ccell*dVdt - Offset;
Sigma = 1.4826*mad(Ic(Quiet),1);
Threshold = Opt.KSigma*Sigma;
Pol = sign(movmean(dVdt, ns(1e-6)));

Above = abs(Ic) > Threshold & ~Masked;
idx = find(Above);
[Start,Stop,Peak,EvPol,Ipk,Q,Base,Width] = deal(zeros(0,1));
if ~isempty(idx)
    br = [1 find(diff(idx) > ns(Opt.GapTime))+1];
    A = idx(br); B = idx([br(2:end)-1 numel(idx)]);
    Pk = arrayfun(@(a,b) max(abs(Ic(a:b))),A,B);
    e = 2;
    while e <= numel(A)
        if (A(e) - B(e-1))*dt < Opt.MergeTime && Pk(e) < Opt.MergeFrac*Pk(e-1)
            B(e-1) = B(e); A(e) = []; B(e) = []; Pk(e) = [];
        else
            e = e + 1;
        end
    end
    PrevStop = 0;
    for e = 1:numel(A)
        s = max(1, A(e) - ns(Opt.PreTime));
        en = min(N, B(e) + ns(Opt.PostTime));
        if e < numel(A), en = min(en, A(e+1) - ns(15e-9)); end
        bw = max([1, PrevStop+1, s - ns(Opt.BaseTime)]):s-1;
        if numel(bw) < 10, bw = max(1, s - ns(Opt.BaseTime)):s-1; end
        b = mean(Ic(bw));
        y = Ic(s:en) - b;
        p = Pol(A(e)); if p == 0, p = sign(sum(y)); end
        [pk,ip] = max(p*y);
        l = ip; while l > 1 && p*y(l) > pk/2, l = l-1; end
        r = ip; while r < numel(y) && p*y(r) > pk/2, r = r+1; end
        Start(end+1,1) = s; Stop(end+1,1) = en; Peak(end+1,1) = s-1+ip; %#ok<AGROW>
        EvPol(end+1,1) = p; Ipk(end+1,1) = p*pk; Q(end+1,1) = trapz(y)*dt; %#ok<AGROW>
        Base(end+1,1) = b; Width(end+1,1) = (r-l)*dt; %#ok<AGROW>
        PrevStop = en;
    end
end
Ev.Ic = Ic; Ev.Offset = Offset; Ev.Sigma = Sigma; Ev.Threshold = Threshold;
Ev.tSwitch = t(iSw); Ev.Masked = Masked; Ev.Pol = Pol; Ev.dVdt = dVdt;
Ev.Events = table(Start, Stop, Peak, (t(Peak)' - t(iSw))*1e6, V(Start)'/1e3, EvPol, Ipk, Q*1e9, Base, Width*1e9, ...
    'VariableNames', {'Start','Stop','Peak','tPeak_us','V_kV','Pol','Ipk_A','Q_nC','Base_A','Width_ns'});
end
