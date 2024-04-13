function [ H, h_true, LAMBDA, partial_LAMBDA ] = partial_LAMBDA_UCyA_f( gamma, AOA, ZOA, UCyA_elements_nor, Nt, Nr_UCA, Nr_ULA, d_H, FL )

    [H, h_true] = gen_chan_specular_UCyA( gamma, AOA, ZOA, UCyA_elements_nor, Nt, Nr_UCA, Nr_ULA, d_H );

    %% LAMBDA
    LAMBDA  = [];
    for jj = 1 : Nt
        lambda_j =[];
        for r = 1 : Nr_UCA * Nr_ULA
            h_rj       = transpose(H(r,jj));
            lambda_rj  = diag(FL*h_rj);
            lambda_j   = [lambda_j; lambda_rj];
        end
        LAMBDA = [LAMBDA lambda_j];
    end


    %% Partial derivative of LAMBDA w.r.t. h_i
    partial_LAMBDA_ll = [];
    for jj = 1 : Nt
        lambda_jj =[];
        for r = 1 : Nr_UCA * Nr_ULA
            lambda_rj_ll = diag(FL(:, 1));
            lambda_jj    = [lambda_jj; lambda_rj_ll];
        end
        partial_LAMBDA_ll = [partial_LAMBDA_ll lambda_jj];
    end
    partial_LAMBDA = partial_LAMBDA_ll;

end