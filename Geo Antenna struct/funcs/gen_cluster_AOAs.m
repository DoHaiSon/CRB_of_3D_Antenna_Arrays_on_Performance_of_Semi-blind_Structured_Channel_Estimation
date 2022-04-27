function cluster_AOAs = gen_cluster_AOAs(ASA, Pn, K, AOA_LOS, N)
% ASA: AOA spread
% Pn: Normalized powers of clusters
% K : Ricean factor
% AOA_LOS: AOA of LOS
% N : number of clusters
% Return: AOAs of clusters
    list = [4, 5, 8, 10, 11, 12, 14, 15, 16, 19, 20];
    if ~ismember(N, list)
        error('Number of clusters should be in this list: 4, 5, 8, 10, 11, 12, 14, 15, 16, 19, 20.');
    end
    
    switch N
        case 4
            C_NLOS_AOA = 0.779;
        case 5
            C_NLOS_AOA = 0.860;
        case 8
            C_NLOS_AOA = 1.018;
        case 10 
            C_NLOS_AOA = 1.090;
        case 11
            C_NLOS_AOA = 1.123;
        case 12
            C_NLOS_AOA = 1.146;
        case 14
            C_NLOS_AOA = 1.190;
        case 15
            C_NLOS_AOA = 1.211;
        case 16
            C_NLOS_AOA = 1.226;
        case 19
            C_NLOS_AOA = 1.273;
        case 20
            C_NLOS_AOA = 1.289;
    end
    
    max_Pn = max(Pn);
    
    C_AOA  = C_NLOS_AOA * (1.1035 - 0.028*K - 0.002*K^2 + 0.0001*K^3);   % for LOS
    
    for i = 1:N
        ln_x  = log10(Pn(i) / max_Pn) / log10(exp(1));
        AOA_1 = 2*( (ASA/1.4)*sqrt(-ln_x) ) / C_AOA;
    end
    
    X_n = unifrnd(-1, 1, N, 1);
    Y_n = normrnd(0, (ASA/7)^2, N, 1);
    
    cluster_AOAs = X_n .* AOA_1 + Y_n + AOA_LOS;
end