function [corr,C] = split_hist_edges(T_pre,T_post,ta,tb,bin_corr,edges)

% calculates the split cross-correlograms of the two neurons.
% T_pre: spike times of the presynaptic neuron
% T_post: spike times of the postsynaptic neuron
% ta: time(seconds) where the correlation window starts (ex: -.05)
% tb: time (seconds) where the window ends (ex: .15)
% bin_corr:number of bins for our correlation window

dt=(tb-ta)/bin_corr;
isi = diff(T_pre);
corr = zeros(length(edges)-1,bin_corr);

for i=1:(length(edges)-1)
    Tlo = edges(i);
    Thi = edges(i+1);
    T = T_pre(find(isi>Tlo & isi<Thi)+1);

    if ~isempty(T)
        [corr(i,:)] = corrFast(T, T_post,ta,tb,bin_corr);
    end    
    C(i)=length(T);
end