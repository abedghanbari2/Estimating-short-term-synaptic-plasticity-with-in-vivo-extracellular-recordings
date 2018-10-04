
function [f,df] = negLogLik_noplas(fit_par,X,Yy,mYy,PSPy)

Np = size(X,2);
bta = fit_par(1:Np);
Mu=(X*bta)';
A = fit_par(Np+1);
gam = PSPy*ones(1,size(X,1))*A;
lam = sigmoidfxn(Mu+gam).*mYy;
lam(lam==1)=1-10e-12;
lam(lam==0)=10e-12;
f = sum(sum(-log(lam).*Yy-log(1-lam).*(1-Yy)));
dbta = sum((lam-Yy)*X,1);
dA = sum(sum( (lam-Yy).*gam/A ));

df = [dbta'; dA];

if any(~isfinite(df))
    f = Inf;
end


