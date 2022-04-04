function [H,h_vec] = gen_chan_specular_ULA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor,Nr,L,N_t)
    % fading, delay, DOA of size (M,Nt)
    H = zeros(Nr,L,N_t);
    M = size(DOA_Theta,1);  

    % Suppose that d/lambda = 1/2
    
    h_vec = [];
    r=0;
    
    for Nr_index=1:Nr
        r=r+1;
        for jj = 1 : N_t
            for ll = 1 : L
                h = 0;
                for mm = 1 : M 
                    r_x = sin(DOA_Theta(mm,jj)) * cos(DOA_Phi(mm,jj));
                    r_y = sin(DOA_Theta(mm,jj)) * sin(DOA_Phi(mm,jj));
                    r_z = cos(DOA_Theta(mm,jj));
                    h = h + fading(mm,jj) ...
                        * sinc((ll-1)-delay(mm,jj)) * exp(-1i*2*pi* ...
                        (position_elements_nor(1, 1, Nr_index)*r_x ...
                        + position_elements_nor(2, 1, Nr_index)*r_y ...
                        + position_elements_nor(3, 1, Nr_index)*r_z));
                end
                H(r,ll,jj) = h;
                h_vec      = [h_vec; h];
            end
        end
    end
end