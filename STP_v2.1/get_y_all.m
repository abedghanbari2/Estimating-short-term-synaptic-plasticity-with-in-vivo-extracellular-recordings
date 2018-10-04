
function [Yy,w,tc,reset] = get_y_all(Tlist,syn,deltat,hyper_params)
% get weights (under the alpha-function) of the post-synaptic spikes
w = zeros(size(Tlist{1}))*NaN;
td = zeros(size(Tlist{1}))*NaN;
q=hyper_params.window_size;
tt = linspace(hyper_params.tblock_min,hyper_params.tblock_max,q+1)';
N=length(Tlist{1});

% restrict to range
deltatp=deltat;
deltatp=deltatp(deltatp(:,1)>hyper_params.tblock_min & deltatp(:,1)< hyper_params.tblock_max,:);
un = unique(deltatp(:,2));
% keyboard
Yy = zeros(q,N);
for i=1:length(un)
    td(un(i)) = min(deltatp(deltatp(:,2)==un(i),1));
    w(un(i)) = syn(td(un(i)));
    [~,ind_tmp] = histc(deltatp(deltatp(:,2)==un(i),1),tt);
    Yy(ind_tmp,un(i)) = 1;
end
w(~isfinite(w))=NaN;
% td(~isfinite(td))=NaN;

% reset
reset = ones(1,N);
for i=2:N
    if Tlist{2}(find(Tlist{2}<Tlist{1}(i),1,'last'))>Tlist{1}(i-1)
        reset(i)=0;
    end
end

tc = (tt(2:end)+tt(1:end-1))/2;
% [~,ind]=histc(td,tt);
% indices = sub2ind([q N], ind(ind~=0), find(isfinite(w)));
% Yy = zeros(q,N);Yy(indices) = 1;