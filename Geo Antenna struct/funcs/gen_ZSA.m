function ZSA = gen_ZSA(fc)
% fc: carrier freq
% Return: ZSA spread in degree

    mean_log_ZSA = -0.1 * log10(1+ fc) + 0.73;
    std_log_ZSA  = -0.04* log10(1+ fc) + 0.34;
    std_log_ZSA  = max(0, std_log_ZSA);
    
    log_ZSA      = normrnd(mean_log_ZSA, std_log_ZSA);
    
    ZSA          = 10.^(log_ZSA);
end