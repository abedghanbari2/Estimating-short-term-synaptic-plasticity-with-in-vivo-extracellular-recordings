function [syn,deltat,estParam,hyper_params]=synapse_xcorr(Tlist,hyper_params)
if isfield(hyper_params,'bin_width')
    bin_width = hyper_params.bin_width;
else
    bin_width = .00002;
end

t = linspace(-hyper_params.coupling_timescale,hyper_params.coupling_timescale,floor(hyper_params.coupling_timescale/bin_width));
[ccgram,deltat] = corrFast(Tlist{1},Tlist{2},min(t),max(t),length(t));
ccgram = ccgram(1:end-1);

XX = getCubicBSplineBasis(linspace(0,1,length(t)-1),hyper_params.baseline_nsplines,0);

% fit model
options=[];
options.method = 'cg';
options.MaxIter = 200;
options.Display = 'off';
f=Inf;
for rr=1:100
    rlat = log(gamrnd(2,.001)+.0001);
    if mod(rr,2)==0
        b01 = [log(nanmean(ccgram)); zeros(size(XX,2)-1,1); rlat; rlat+randn(1); randn(1)];
    else
        b01 = [log(nanmean(ccgram)); randn(size(XX,2)-1,1)/5; rlat; rlat+randn(1); randn(1)];
    end
    [brr,frr] = minFunc(@lossGlmAlpha,b01,options,XX,ccgram,t(1:end-1)');
    if frr<f
        b1=brr;
        f=frr;
    end
end

options.MaxIter=2000;
[b1,~] = minFunc(@lossGlmAlpha,b1,options,XX,ccgram,t(1:end-1)');

estParam.syn_params = b1(end-2:end);
estParam.syn_params(1:2) = exp(estParam.syn_params(1:2));
syn = @(ts) max(0,ts-estParam.syn_params(1))/estParam.syn_params(2).*exp(1-max(0,ts-estParam.syn_params(1))/estParam.syn_params(2));
estParam.b1=b1;
% basic window
thr = hyper_params.thr;
% hyper_params.tblock_min = min(t(abs(syn(t))>thr*abs(max(syn(t)))));
hyper_params.tblock_min = 0;
hyper_params.tblock_max = max(t(abs(syn(t))>thr*abs(max(syn(t)))));
