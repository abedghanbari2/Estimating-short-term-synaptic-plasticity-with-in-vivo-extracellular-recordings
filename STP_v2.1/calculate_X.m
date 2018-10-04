function [X,history,hyper_params] = calculate_X(Tlist,hyper_params)

% spline basis over time
mint = min(cellfun(@min,Tlist));
maxt = max(cellfun(@max,Tlist));
tt = (Tlist{1}-mint)/(maxt-mint);
if ~isfield(hyper_params,'nonstationary_nsplines')
    hyper_params.nonstationary_nsplines = ceil((maxt-mint)/hyper_params.XTimeSpan);
end

Xfr = getCubicBSplineBasis(tt,hyper_params.nonstationary_nsplines,0);
Xfr = Xfr(:,2:end);

% history basis
history.nbas = hyper_params.hist_nbs;
history.maxt = hyper_params.hist_mxt;
history.pow = 2;
dk = 1/(history.nbas-1);
history.knots = linspace(-1/4*dk-2*dk,1,history.nbas+4);
history.knots = (abs(history.knots)).^history.pow.*sign(history.knots);
history.knots = history.knots/max(history.knots)*history.maxt;
history.spline = fastBSpline(history.knots,ones(history.nbas,1));
N = size(tt,1);
XH = zeros(N,history.nbas);
for i=1:N
    k=find(Tlist{2}<Tlist{1}(i) & Tlist{2}>Tlist{1}(i)-history.maxt,1,'last');
    if ~isempty(k)
        XH(i,:) = (history.spline.getBasis(Tlist{1}(i)-Tlist{2}(k)));
    end
end
XH=XH(:,2:end);
X = [Xfr XH];
% X = [ XH];


