
function [f,df,lam] = lossGlmAlpha(b,X,y,t)

p = b((size(X,2)+1):end);
p(1:2)=exp(p(1:2));
b = b(1:size(X,2));

t(t<p(1))=p(1);
syn = p(3)*(t-p(1))/p(2).*exp(1-(t-p(1))/p(2));

lam = exp(X*b + syn);
f = -nansum(y.*log(lam+(lam==0))-lam);

fb = sum(b(2:end).^2);
f = f+fb;

lam_err = lam-y;
lam_err(~isfinite(lam_err))=0;
db = X'*lam_err;
dp = [nansum((syn/p(2)-syn./(t-p(1))).*lam_err);...
    ((-syn/p(2)+syn.*(t-p(1))/p(2).^2))'*lam_err;...
    syn'*lam_err/p(3)];
dp(1:2)=dp(1:2) .* p(1:2);
db = db+2*[0; b(2:end)];
df = [db; dp];
