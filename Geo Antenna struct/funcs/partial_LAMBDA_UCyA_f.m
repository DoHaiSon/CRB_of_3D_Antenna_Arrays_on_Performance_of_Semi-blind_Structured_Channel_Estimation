function [ H, h_true, LAMBDA, partial_LAMBDA ] = partial_LAMBDA_UCyA_f( fading, delay, DOA_Phi, DOA_Theta, position_elements_nor, Nr_UCA, Nr_ULA, L, Nt, FL )

    [H, h_true] = gen_chan_specular_UCyA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor,Nr_UCA,Nr_ULA,L,Nt);

    %% LAMBDA
    LAMBDA  = [];
    for jj = 1 : Nt
        lambda_j =[];
        for r = 1 : Nr_UCA * Nr_ULA
            h_rj       = transpose(H(r,:,jj));
            lambda_rj  = diag(FL*h_rj);
            lambda_j   = [lambda_j; lambda_rj];
        end
        LAMBDA = [LAMBDA lambda_j];
    end


    %% Partial derivative of LAMBDA w.r.t. h_i
    partial_LAMBDA  = cell(1, L);
    for ll = 1 : L
        partial_LAMBDA_ll = [];
        for jj = 1 : Nt
            lambda_jj =[];
            for r = 1 : Nr_UCA * Nr_ULA
                lambda_rj_ll = diag(FL(:,ll));
                lambda_jj    = [lambda_jj; lambda_rj_ll];
            end
            partial_LAMBDA_ll = [partial_LAMBDA_ll lambda_jj];
        end
        partial_LAMBDA{1,ll} = partial_LAMBDA_ll;
    end

end