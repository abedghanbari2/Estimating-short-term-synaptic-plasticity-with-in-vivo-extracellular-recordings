
function [f,df] = negLogLik(fit_par,ts,X,Yy,mYy,PSPy,reset,model_str,printresults)
% keyboard
% td : time difference between first post synaptic spike after a
% presynaptic spike in the window

% reset : {1: if there was a postsynaptic spike after previous presynaptic
% spike}  { 0: if no postsynaptic spike }

etmparams = fit_par(end-4:end);
Np = size(X,2);
bta = fit_par(1:Np);
Mu=(X*bta)';
A = fit_par(Np+1);
etmparams_n( [1 2 5]) = exp(etmparams( [1 2 5]));
etmparams_n( [3 4]) = 1./(1+exp(-etmparams( [3 4])));

p = eTM_v2(ts,reset,etmparams_n,model_str);
p = p/mean(p);
gam = PSPy*p'*A;
lam = sigmoidfxn(Mu+gam).*mYy;
% lam(lam==1)=1-10e-12;
% lam(lam==0)=10e-12;

% f = sum(sum(-(Mu+gam).*Yy-log(1-lam))) + sum(etmparams.^2);
f = sum(sum(-log(lam+(lam==0)).*Yy-log(1-lam + (lam==1)).*(1-Yy))) + sum(etmparams.^2);

if nargout>1
    dbta = sum((lam-Yy)*X,1);
    dA = sum(sum( (lam-Yy).*gam/A ));

    % DFUfS
    e = 10e-8;
    detm=zeros(size(etmparams));

    for i=1:length(etmparams)
        detmparams = zeros(size(etmparams));
        detmparams(i) = e;
        fit_par1 = fit_par;
        fit_par1(end-4:end) = etmparams + detmparams;
        f1 = negLogLik(fit_par1,ts,X,Yy,mYy,PSPy,reset,model_str,printresults);

        fit_par2 = fit_par;
        fit_par2(end-4:end) = etmparams - detmparams;
        f2 = negLogLik(fit_par2,ts,X,Yy,mYy,PSPy,reset,model_str,printresults);

        detm(i) = -(f2-f1)/2/e;
    end

    df = [dbta'; dA; detm];

    if printresults
        for i = 1:5
            fprintf('%4.2f \t',etmparams_n(i))
        end
        fprintf('%4.2f \n',f)
    end

    if any(~isfinite(df))
        f = Inf;
    end

end
