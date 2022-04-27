function cluster_ZOAs = gen_cluster_ZOAs(ZSA, Pn, K, ZOA_LOS, N)
% ZSA: ZOA spread
% Pn: Normalized powers of clusters
% K : Ricean factor
% ZOA_LOS: ZOA of LOS
% N : number of clusters
% Return: ZOAs of clusters
    list = [8, 10, 11, 12, 15, 19, 20];
    if ~ismember(N, list)
        error('Number of clusters should be in this list: 8, 10, 11, 12, 15, 19, 20.');
    end
    
    switch N
        case 8
            C_NLOS_ZOA = 0.889;
        case 10 
            C_NLOS_ZOA = 0.957;
        case 11
            C_NLOS_ZOA = 1.031;
        case 12
            C_NLOS_ZOA = 1.104;
        case 15
            C_NLOS_ZOA = 1.1088;
        case 19
            C_NLOS_ZOA = 1.184;
        case 20
            C_NLOS_ZOA = 1.178;
    end
    
    max_Pn = max(Pn);
    
    C_AOA  = C_NLOS_ZOA * (1.3086 + 0.0339*K - 0.0077*K^2 + 0.0002*K^3);   % for LOS
    
    for i = 1:N
        ln_x  = log10(Pn(i) / max_Pn) / log10(exp(1));
        ZOA_1 = ZSA *ln_x / C_AOA;
    end
    
    X_n = unifrnd(-1, 1, N, 1);
    Y_n = normrnd(0, (ZSA/7)^2, N, 1);
    
    cluster_ZOAs = X_n .* ZOA_1 + Y_n + ZOA_LOS;
end