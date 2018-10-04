function result = estimateSTP(Tlist,models)

% models = {'noSTP','IntegrationOnly','Integration_subThreshold','FacilitationOnly','DepressionOnly','TM','eTM_woIntegration','eTM_subThreshold','eTM_Integration'};

hyper_params.coupling_timescale = 0.005;
hyper_params.baseline_nsplines = 4;
hyper_params.thr =.05;
hyper_params.hist_nbs = 4; % number of history basis
hyper_params.hist_mxt = .01; % sec
hyper_params.window_size = 50;
hyper_params.XTimeSpan = 50; %sec
hyper_params.nrr = 5;

        
for i=1:length(Tlist)
    Tlist{i}=sort(Tlist{i}+(rand(size(Tlist{i}))-.5)*0.0001);
end

% synapse (alpha_func) estimation
[syn,deltat,synParams,hyper_params_tmp]=synapse_xcorr(Tlist,hyper_params);

% output - zero OR one for each presynaptic spike
[Yy,w,tc,reset] = get_y_all(Tlist,syn,deltat,hyper_params_tmp);

% Yy = Yy(:,randperm(size(Yy,2))); % shuffle

% covariates for each presynaptic spike
[X,history,hyper_params_tmp] = calculate_X(Tlist,hyper_params_tmp);
mYy = mask_Yy(Yy);

% estimate params
[stpParams,ll] = fitSTP(Tlist{1},X,Yy,mYy,syn(tc),reset,models,hyper_params_tmp);

for m = 1:length(models)
    [p1, Integ, DF,u,R] = eTM_v2(Tlist{1},reset,stpParams(end-4:end,m),models{m});
    p=p1/mean(p1);
    mu = X*stpParams(1:size(X,2),m);
    nu = ones(hyper_params_tmp.window_size,1)*mu'+stpParams(end-5,m)*syn(tc)*p';
%     nu = mean(ones(hyper_params_tmp.window_size,1)*mu')+stpParams(end-5,m)*syn(tc)*p';
    lam = sigmoidfxn(nu);
    
    antiLam = 1 - lam;
    PSPcm = lam(1,:);
    for i = 2:size(lam,1)
        PSPcm = PSPcm + lam(i,:).*prod(antiLam(1:i-1,:),1);
    end
    PSP(:,m)  = PSPcm' ./ (PSPcm'+prod(antiLam,1)');
%     PSP(:,m)  = sum(lam,1);
    urs_avg(:,m) = mean([Integ u R]);
    m_x(m) =  mean(X*stpParams(1:end-6,m));
    m_xfr(m) =  mean(X(:,1:hyper_params_tmp.nonstationary_nsplines)*stpParams(1:hyper_params_tmp.nonstationary_nsplines,m));

end


% output
result.history = history;
result.synParams = synParams;
result.stpParams = stpParams;
result.ll = ll;
result.models = models;
result.hyper_params = hyper_params_tmp;
result.PSP = PSP;
result.w = w;
result.m_x = m_x;
result.m_xfr = m_xfr;
result.urs_avg = urs_avg;
result.syn = syn;


result.lam = lam; % last lam
result.Yy = Yy;