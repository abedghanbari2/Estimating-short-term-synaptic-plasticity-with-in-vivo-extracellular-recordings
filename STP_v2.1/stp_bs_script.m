% clc;clear variables;close all;
restoredefaultpath;

addpath(genpath('STP_v2.1'))
addpath(genpath('data_remote'))

models = {'noSTP','IntegrationOnly','Integration_subThreshold','FacilitationOnly','DepressionOnly','TM','eTM_woIntegration','eTM_subThreshold','eTM_Integration'};

% load data
nrep = (feature('numCores'));
parpool(nrep)

nfold = 12;
parfor irep =  1:nrep
    Tlist = load_data_invivo(ni);
    T_chunk = floor(max(max(Tlist{1}),max(Tlist{2}))/nfold);
    Tlist_cv = bootstrap_population(Tlist,T_chunk,nfold);
    result{irep} = estimateSTP(Tlist_cv,models)
end
save(['result',num2str(ni),'.mat'])
