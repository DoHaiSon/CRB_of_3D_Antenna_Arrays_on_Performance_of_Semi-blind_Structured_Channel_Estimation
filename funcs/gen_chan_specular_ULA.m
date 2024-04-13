function [H, h_vec] = gen_chan_specular_ULA(gamma, AOA, ZOA, ULA_elements_nor,  Nt, Nr, d_H)
    % fading, delay, DOA of size (M,Nt)
    H = zeros(Nr, Nt);

    h_vec = [];    
    for Nr_index=1:Nr
        for jj = 1 : Nt
            h = 0;
            for nn = 1 : d_H 
                r_x = sin(ZOA(nn,jj)) * cos(AOA(nn,jj));
                r_y = sin(ZOA(nn,jj)) * sin(AOA(nn,jj));
                r_z = cos(ZOA(nn,jj));
                h = h + gamma(nn,jj) ...
                    * exp(-1i*2*pi* ...
                     (ULA_elements_nor(1, 1, Nr_index)*r_x ...
                    + ULA_elements_nor(2, 1, Nr_index)*r_y ...
                    + ULA_elements_nor(3, 1, Nr_index)*r_z));
            end
            H(Nr_index, jj) = h;
            h_vec = [h_vec; h];
        end
    end
end