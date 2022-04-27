function cluster_delays = gen_cluster_delays(rt, DS, N)
% rt: Delay scaling param
% DS: Delay spread
% N : number of clusters
% K : Ricean factor
% Return: delays of clusters

    t_n_1 = -rt * DS .* (log10(rand(N, 1) / log10(exp(1))));
    
    t_n = sort(t_n_1 - min(t_n_1));
    
    cluster_delays = t_n;
end