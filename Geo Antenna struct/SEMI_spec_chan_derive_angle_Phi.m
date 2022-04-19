function Br_angle_Phi= SEMI_spec_chan_derive_angle_Phi(mag,phase,delay,DOA_Phi,DOA_Theta,position_elements_nor,L,N,M,Nt)
%Nt   % number of transmit antennas
%Nr   % number of receive antennas
%L    % channel order
%N    % number of TCs
%M    % number of sub-paths per TC

    Br_angle_Phi_tmp = zeros(N,M,L,Nt);
    for jj = 1 : Nt
        for nn = 1 : N
            for mm = 1 : M
                for l = 1 : L
                    r_x = sin(DOA_Theta(nn,mm,jj)) * cos(DOA_Phi(nn,mm,jj));
                    r_y = sin(DOA_Theta(nn,mm,jj)) * sin(DOA_Phi(nn,mm,jj));
                    r_z = cos(DOA_Theta(nn,mm,jj));
                    Br_angle_Phi_tmp(nn,mm,l,jj)=mag(nn,mm,jj)* exp(1i * phase(nn,mm,jj)) * sinc((l-1)-delay(nn,mm,jj)) ...
                        * exp(-1i*2*pi*(position_elements_nor(1)*r_x + position_elements_nor(2)*r_y + position_elements_nor(3)*r_z)) ...
                        * (-1i*2*pi)*(-position_elements_nor(1) * sin(DOA_Theta(mm,jj)) * sin(DOA_Phi(mm,jj)) + ...
                        + position_elements_nor(2) * sin(DOA_Theta(mm,jj)) * cos(DOA_Phi(mm,jj))  + position_elements_nor(3) * r_z);
                end
            end
        end
    end
    Br_angle_Phi_tmp1=cell(1,Nt);
    for jj = 1 : Nt
        Br_angle_Phi_tmp1{1,jj}=Br_angle_Phi_tmp(:,:,:,jj);
    end
    Br_angle_Phi=blkdiag(Br_angle_Phi_tmp1{:});
end

