% clc;
clear variables;
% close all;
restoredefaultpath;

if ispc
    addpath(genpath('C:\Users\abg14006\Google Drive\stp in vivo\code\STP_v2.1'))
    addpath(genpath('C:\Users\abg14006\Google Drive\stp in vivo\data\data_remote'))
else
    addpath(genpath('stpSpike_v2'))
    addpath(genpath('data_remote'))
end

% hyper-params
hyper_params.coupling_timescale = 0.005;
hyper_params.baseline_nsplines = 4;
hyper_params.thr =.2;
hyper_params.hist_nbs = 4; % number of history basis
hyper_params.hist_mxt = .01; % sec
hyper_params.window_size = 25;
hyper_params.nrr = 7;

% load data
figure(1),clf

for neurons = [2:4 6 11]

        Tlist = load_data_invivo(neurons);
        
        for i=1:length(Tlist)
            Tlist{i}=sort(Tlist{i}+(rand(size(Tlist{i}))-.5)*0.0001);
        end

        % synapse (alpha_func) estimation
        [syn,deltat,estParam,hyper_params_tmp]=synapse_xcorr(Tlist,hyper_params);

        % output - zero OR one for each presynaptic spike
        [Yy,w,tc,reset] = get_y_all(Tlist,syn,deltat,hyper_params_tmp);

        % covariates for each presynaptic spike
        [X,history,hyper_params_tmp] = calculate_X(Tlist,hyper_params_tmp);
        mYy = mask_Yy(Yy);

        % estimate params
        %           1           2                       3                   4                   5           6           7                   8               9
        % models = {'noSTP','IntegrationOnly','Integration_subThreshold','FacilitationOnly','DepressionOnly','TM','eTM_woIntegration','eTM_subThreshold','eTM_Integration'};
        models = {'eTM_Integration'};
        [stpParams{neurons},ll] = fitSTP(Tlist{1},X,Yy,mYy,tc,syn,reset,models,hyper_params_tmp);

        % plot
        m = 1;
        ss([2:4 6 11]) = 1:5;
        col_plot = 5;

        [p1, Integ, DF,u,R] = eTM_v2(Tlist{1},reset,stpParams{neurons}(end-4:end,m),models{m});
        p=p1/mean(p1);
        mu = X*stpParams{neurons}(1:size(X,2),m);
        nu = ones(hyper_params_tmp.window_size,1)*mu'+stpParams{neurons}(end-5,m)*syn(tc)*p';

        lam = sigmoidfxn(nu);
        antiLam = 1 - lam;
        PSPcm = lam(1,:);
        for i = 2:size(lam,1)
            PSPcm = PSPcm + lam(i,:).*prod(antiLam(1:i-1,:),1);
        end
        PSPcm  = PSPcm';

        subplot(4,col_plot,ss(neurons))
        hold on
        isi1 = logspace(-3,1,100);
        ppr=[];
        for i=1:length(isi1)
            train = [1:.01:2 10 10+isi1(i)];
            [PSP, Integ, DF]=eTM_v2(train,0*train+1,stpParams{neurons}(end-4:end,m),models{m});
            ppr(i,1)=PSP(end)/PSP(end-1);
            ppr(i,2) = Integ(end)/PSP(end-1);
            ppr(i,3)=DF(end)/PSP(end-1) ;
        end
        plot(log10(isi1),ppr(:,1),'r')
        plot(log10(isi1),ppr(:,2),'g')
        plot(log10(isi1),ppr(:,3),'b')
        ylim([0 1.1*max(ppr(:,1))])
        
        isi = [median(diff(Tlist{1})); diff(Tlist{1})];
        logisi = log10(isi);
        idx=isfinite(logisi);
        logisi=logisi(idx);
        isi_bb = min(isi)-eps; 
        isi_ub = min(max(diff(Tlist{1})),2)+eps;
        edges = logspace(log10(isi_bb),log10(isi_ub),15);  
        xlb = max(log10(isi_bb),min(logisi)); xub=min(log10(isi_ub),max(logisi));
        xlp = linspace(xlb,xub,50);

        tcorr1 = linspace(-.001,.005,100);
        [corr1,C] = split_hist_edges(Tlist{1},Tlist{2},min(tcorr1),max(tcorr1),length(tcorr1),edges);
        dat = sum(corr1(:,tcorr1>hyper_params_tmp.tblock_min & tcorr1<hyper_params_tmp.tblock_max),2)./C';

        mu = (X*stpParams{neurons}(1:size(X,2),m));
        A = stpParams{neurons}(size(X,2)+1,m);
        q=hyper_params_tmp.window_size;
        tt = linspace(hyper_params_tmp.tblock_min,hyper_params_tmp.tblock_max,q+1)';
        tc = (tt(2:end)+tt(1:end-1))/2;
        dt=mean(diff(tt));

        subplot(4,col_plot,col_plot+ss(neurons))

        ylp = linspace(0,max(max(dat),max(PSPcm(idx)))+10e-6,100); % change the max so it won't blow up on y axis
        ns = histcounts2(PSPcm(idx),logisi,ylp,xlp);
        ns = (bsxfun(@rdivide,ns,sum(ns)));
        ym = ns'*(ylp(1:end-1)+mean(diff(ylp))/2)';
        % plot(xlp(1:end-1)+mean(diff(xlp))/2,(ym),'b','LineWidth',2)
        [p1, s] = polyfit(log10(isi(idx)),(PSPcm(idx)),10);
        x = linspace(min(log10(isi(idx))),max(log10(isi(idx))),100);
        [y_fit, delta] = polyval(p1, x, s);
        Xx=[x,fliplr(x)];
        Y=[(y_fit + 2*delta),fliplr((y_fit - 2*delta))];
        h = fill(Xx,Y,'b','edgecolor','none');
        set(h,'facealpha',.1)
        hold on
        plot(x,(y_fit))
        scatter(log10(isi),(PSPcm),'k.','markeredgealpha',.05)
        box off; set(gca,'TickDir','out')
        axis tight
        xlim([log10(.001) log10(3)])
        plot(log10(edges(1:end-1))+mean(diff(log10(edges(1:end-1))))/2,dat,'r-.')
        ylim([0 2.1*max(dat)])

        subplot(4,col_plot,2*col_plot+ss(neurons))
        pp =PSPcm;
        xrg = linspace(min(pp),max(pp),50);
        [ym,~,ns] = bindata(double(isfinite(w)),pp,xrg);
        plot(xrg,ym,'o')
        hold on
        errorbar(xrg,ym,sqrt(ym.*(1-ym)./ns),'.')
        plot(xrg,xrg)

        subplot(4,col_plot,3*col_plot+ss(neurons))
        xrg = linspace(min(logitfxn(pp)),max(logitfxn(pp)),50);
        [ym,~,ns] = bindata(double(isfinite(w)),logitfxn(pp),xrg);

        plot(xrg,ym,'o'); hold on
        errorbar(xrg,ym,sqrt(ym.*(1-ym)./ns),'.')
        plot(xrg,1./(1+exp(-xrg)))
        hold off; axis tight

end
