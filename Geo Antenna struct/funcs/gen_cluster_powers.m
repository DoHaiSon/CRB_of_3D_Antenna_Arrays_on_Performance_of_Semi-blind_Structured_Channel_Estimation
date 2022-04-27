function cluster_powers = gen_cluster_powers(delay_n, rt, DS, gamma, N)
% delay_n: delays of clusters
% rt: Delay scaling param
% DS: Delay spread
% N : number of clusters
% gamma: SF's std
% Return: powers of clusters

    Z_n   = normrnd(0, gamma, N, 1);
    P_n_1 = exp(-delay_n .* ((rt - 1) / (rt * DS))) .* 10.^(-Z_n / 10);
    
    Pr   = sum(P_n_1);
    
    cluster_powers = P_n_1 / Pr;
end