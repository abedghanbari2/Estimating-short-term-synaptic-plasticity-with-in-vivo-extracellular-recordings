function [population,dt,T,delay,coupling] = load_data_invivo(data_category,varargin)

switch data_category
    case 1
        load('July 14 A.mat')
        population{1} = LGN;
        population{2} = CTX;
        dt = .001;
        delay = [5 100 200]; % [coupling history common_input] ms
        coupling.nfilt = 6;
        T = 2200;

    case 2
        load('July13A.mat')
        population{1} = LGN;
        population{2} = CTX;
        dt = .001;
        delay = [5 100 200]; % [coupling history common_input] ms
        coupling.nfilt = 6;
        T = 2900;
    case 3
        load('Tuning_G1003C82R14.mat')
        population{1} = Tuning.Raw.Continuous.STsEPSP/Tuning.SamplingRate;
        population{2} = Tuning.Raw.Continuous.STs(Tuning.Raw.Continuous.isAP)/Tuning.SamplingRate; 
        dt = .0001;
        delay = [3 100 100]; % [coupling history common_input] ms
        coupling.nfilt = 10;
        T = 350;
    case 4
        load('crcns-ch-epfl-2009-challengeD.mat')
%         rng(1);
        population{1} = Tlist{1};
        population{2}=Tlist{2};
        for i=1:length(population{2})
            if ~isempty(find(population{1}==population{2}(i), 1))
                population{2}(i) = population{2}(i)+.001;
            end
        end
     
        population{2} = sort(population{2}+(randn(size(population{2})))*0.00005); 
        dt = .001;

    case 5
        load('gor01-6-13_multisess.mat')
        population{1} = Tlist{35};
        population{2} = Tlist{18}+.002;
        dt = .001;
    case 6
        load('Tlist_abed.mat')
        population{1} = Tlist{1,800};
        population{2} = Tlist{1,789};
        dt = .001;
        
    case 7
        load('Tlist_abed.mat')
%         population{1} = Tlist{2,654};
%         population{2} = Tlist{2,689};
        population{1} = Tlist{1,167};
        population{2} = Tlist{1,106};
        dt = .001;
        
    case 8
        load('Spont_Tuning_NatSounds_G1007C7_v2.mat')
        SamplingRate = 97.65625e3;
        population{1} = NatSounds.continuous.STsEPSP/SamplingRate;
        population{2} = NatSounds.continuous.STs(NatSounds.continuous.isAP)/SamplingRate; 
        dt = .0001;
        
    case 9
        load('Spont_Tuning_NatSounds_G1007C7_v2.mat')
        SamplingRate = 97.65625e3;
        population{1} = Spont.STsEPSP/SamplingRate;
        population{2} = Spont.STs(Spont.isAP)/SamplingRate; 

    case 10
        load('Spont_Tuning_NatSounds_G1007C7_v2.mat')
        SamplingRate = 97.65625e3;
        population{1} = Tuning.Continuous.STsEPSP/SamplingRate;
        population{2} = Tuning.Continuous.STs(Tuning.Continuous.isAP)/SamplingRate; 
    case 11
        load('Mar23c1,2,4_Herc_spikes.mat')
        population{1} = Tlist{1} - floor(min([Tlist{1}; Tlist{2}]));
        population{2} = Tlist{2} - floor(min([Tlist{1}; Tlist{2}]));
        delay = [5 100 200]; % [coupling history common_input] ms
        coupling.nfilt = 6;

    otherwise
        disp('choose between 1 2 3')
end

