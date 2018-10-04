function [PSP, Integ, DF,utmp,Rtmp] =eTM_v2(spk,reset,params,model) %#codegen

D = params(1);
F = params(2);
U = params(3);
f = params(4);
S = params(5);
switch model 
    case 'noSTP'
        D = 0;
        F = 0;
        U = 1;
        f = U;
        reset = 0*reset;
    case 'IntegrationOnly'
        D = 0;
        F = 0;
        U = 1;
        f = U;
    case 'Integration_subThreshold'
        D = 0;
        F = 0;
        U = 1;
        f = U;
        reset = ones(size(reset));
    case 'FacilitationOnly'
        D = 0;
        f = U;
    case 'DepressionOnly'
        F=0;
        f = U;
    case 'TM'
        f = U;
    case 'eTM_woIntegration'
        reset = 0*reset;
    case 'eTM_subThreshold'
        reset = ones(size(reset));
    case 'eTM_Integration'
        %
    case 'singleISI'
        f = U;
        reset = zeros(size(reset));

    otherwise
        error('specify model')
end

L = length(spk);
R = zeros(L,1);
u = zeros(L,1);
Rtmp = zeros(L,1);
utmp = zeros(L,1);
Integ = zeros(L,1);
PSP = zeros(L,1);
DF = zeros(L,1);

R(1) = 1;
u(1) = U;

disi = exp(-diff(spk)/D);
fisi = exp(-diff(spk)/F);
sisi = exp(-diff(spk)/S);

PSP(1) = R(1)*u(1);
for i = 1 : L-1

    if strcmp(model,'singleISI')
        R(i+1) = 1 - U * disi(i);
        u(i+1) =  U + U * (1 - U) *  fisi(i);
        Integ(i+1) = U * sisi(i);

    else
    
        R(i+1) = 1 - ( 1 - R(i) * (1-u(i)) ) * disi(i);
        u(i+1) =  U + ( u(i) + f * ( 1 - u(i)) - U  )*  fisi(i);
        Integ(i+1) = reset(i+1)* PSP(i) * sisi(i);
    
    end
    
    Rtmp(i+1)=R(i+1);
    utmp(i+1)=u(i+1);
    
    if isfinite(R(i+1)) 
        PSP(i+1) =  Integ(i+1) + R(i+1) * u(i+1);
%         PSP(i+1) =  R(i+1) * u(i+1);
    else
        PSP(i+1) =  Integ(i+1);
    end
    DF(i+1) = R(i+1) * u(i+1);
    if ~isfinite(PSP(i+1))
        keyboard
    end
end
PSP(1) = median(PSP(2:end));

















