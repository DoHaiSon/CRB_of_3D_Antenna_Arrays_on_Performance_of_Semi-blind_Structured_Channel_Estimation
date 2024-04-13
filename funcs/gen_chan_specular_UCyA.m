function [H,h_vec] = gen_chan_specular_UCyA( gamma, AOA, ZOA, UCyA_elements_nor, Nt, Nr_UCA, Nr_ULA, d_H )
    % fading, delay, DOA of size (M,Nt)
    H = zeros(Nr_ULA*Nr_UCA, Nt);
    
    h_vec = [];
    r=0;
    for Nr_ULA_index=1:Nr_ULA   
        for Nr_UCA_index=1:Nr_UCA
            r=r+1;
            for jj = 1 : Nt
                h = 0;
                for nn = 1 : d_H 
                    r_x = sin(ZOA(nn,jj)) * cos(AOA(nn,jj));
                    r_y = sin(ZOA(nn,jj)) * sin(AOA(nn,jj));
                    r_z = cos(ZOA(nn,jj));
                    h = h + gamma(nn,jj) ...
                        * exp(-1i*2*pi* ...
                         (UCyA_elements_nor(1, Nr_ULA_index, Nr_UCA_index)*r_x ...
                        + UCyA_elements_nor(2, Nr_ULA_index, Nr_UCA_index)*r_y ...
                        + UCyA_elements_nor(3, Nr_ULA_index, Nr_UCA_index)*r_z));
                end
                H(r, jj) = h;
                h_vec    = [h_vec; h];
            end
        end
    end
end
