function [stpParamsAll,ll] = fitSTP(ts,X,Yy,mYy,syn_tc,reset,models,hyper_params)


options_rr=[];
options_rr.method = 'lbfgs';
% options_rr.MaxIter = 500;
nrr=hyper_params.nrr;
[b00,~] = minFunc(@negLogLik_noplas,randn(size(X,2)+1,1),options_rr,X,Yy,mYy,syn_tc);
% checkgrad('negLogLik_noplas',randn(size(X,2)+1,1),10e-12,X,Yy,mYy,syn_tc)
% keyboard
options_rr.MaxIter = 25;
options_rr.Display = 'off';

for m = 1:length(models)
    f=Inf;
    fprintf(['\n \n',models{m}, '\n'])
    for rr=1:nrr
        b_rr = [b00'+randn(size(b00'))/2 , [-1 -1 0 0 -3] + 5*randn(1,5)]' ;
        [synParams_tmp,fval] = minFunc(@negLogLik,b_rr,options_rr,ts,X,Yy,mYy,syn_tc,reset,models{m},0);
    %     checkgrad('negLogLik',b_rr,10e-12,ts,X,Yy,mYy,syn_tc,reset,models{m},0)
    %     keyboard
        if fval<f
            b0=synParams_tmp';
            f=fval;
            fprintf('repeat %d / %d fvalmin %f \n',rr,nrr,fval)
        end
    end

    options=[];
    options.method = 'lbfgs';
    % options.optTol = 1e-2;
    % options.progTol = 1e-4;
    options.MaxIter = 1000;
    options.MaxFunEvals = 2000;
    options.Display = 'off';
%     fprintf('D \t\t F \t\t U \t\t f \t\t S \t\t fval \n')
    [synParamsT,fval] = minFunc(@negLogLik,b0',options,ts,X,Yy,mYy,syn_tc,reset,models{m},1);
    stpParams=synParamsT;
    stpParams( [1 2 5] +size(X,2)+1) = exp(stpParams( [1 2 5] +size(X,2)+1));
    stpParams( [3 4] +size(X,2)+1) = 1./(1+exp(-stpParams( [3 4] +size(X,2)+1)));

    stpParamsAll(:,m) = stpParams;
    ll(m) = -fval; % min neg log-lik
end
