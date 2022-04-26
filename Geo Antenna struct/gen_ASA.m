function ASA = gen_ASA(fc)
% fc: carrier freq
% Return: ASA spread in degree

    mean_log_ASA = -0.08 * log10(1+ fc) + 1.73;
    std_log_ASA  = 0.014 * log10(1+ fc) + 0.28;
    
    log_ASA      = normrnd(mean_log_ASA, std_log_ASA);
    
    ASA          = 10.^(log_ASA);
end