function DS = gen_DS(fc)
% fc: carrier freq
% Return: delay spread

    mean_log_DS = -0.24 * log10(1+ fc) - 7.14;
    std_log_DS  = 0.34;
    
    log_DS      = normrnd(mean_log_DS, std_log_DS);
    
    DS          = 10.^(log_DS);
end