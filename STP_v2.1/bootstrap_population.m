function [T_bs] = bootstrap_population(population,T_chunk,n_bs_chunks)
%
% T_chunk = 200;
% n_bs_chunks = 5;
%


T = max( max(population{1}) , max(population{2}) );

num_chunks = floor(T/T_chunk);

for i_chunks = 1:num_chunks
    for i_neuron = 1:2
        population_chunk{i_chunks}{i_neuron} = ...
            population{i_neuron}( population{i_neuron} < T_chunk*i_chunks...
            & T_chunk*(i_chunks-1) < population{i_neuron} ) - T_chunk*(i_chunks-1);
    end
end

bs_nlist = randi(num_chunks,1,n_bs_chunks);

for i_neuron=1:2;T_bs{i_neuron}=[];end
for i_bs_chunks = 1:n_bs_chunks
    for i_neuron = 1:2
        T_bs{i_neuron} = [T_bs{i_neuron};...
            (i_bs_chunks-1)*T_chunk + population_chunk{bs_nlist(i_bs_chunks)}{i_neuron}];
    end
end
