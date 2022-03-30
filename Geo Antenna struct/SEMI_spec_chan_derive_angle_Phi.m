function Br_angle_Phi= SEMI_spec_chan_derive_angle_Phi(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor,L,M,Nt)

%Nt = 4;    % number of transmit antennas
%Nr = 4;    % number of receive antennas
%L   = 4;    % channel order
%M   = 2;    % Number of multipaths (assumption: M  = L)   
%fading = rand(M,Nt)+1i*rand(M,Nt);
%fading = rand(M,Nt);
%delay  = rand(M,Nt);
%DOA    = pi * rand(M,Nt);
%d_nor=1/2;
%Nr_index=1;
%Br_fading=zeros()
    Br_angle_Phi_tmp = zeros(M,L,Nt);
    for jj = 1 : Nt
        for mm = 1 : M
            for l = 1 : L
                r_x = sin(DOA_Theta(mm,jj)) * cos(DOA_Phi(mm,jj));
                r_y = sin(DOA_Theta(mm,jj)) * sin(DOA_Phi(mm,jj));
                r_z = cos(DOA_Theta(mm,jj));
                Br_angle_Phi_tmp(mm,l,jj)=fading(mm,jj) * sinc((l-1)-delay(mm,jj)) ...
                    * exp(-1i*2*pi*(position_elements_nor(1)*r_x + position_elements_nor(2)*r_y + position_elements_nor(3)*r_z)) ...
                    * (-1i*2*pi)*(-position_elements_nor(1) * sin(DOA_Theta(mm,jj)) * sin(DOA_Phi(mm,jj)) + ...
                    + position_elements_nor(2) * sin(DOA_Theta(mm,jj)) * cos(DOA_Phi(mm,jj))  + position_elements_nor(3) * r_z);
                %             Br_angle_Phi_tmp(mm,l,jj)=fading(mm,jj)*exp(-1i*2*pi*R_nor*sin(DOA_Theta(mm,jj))*cos(DOA_Phi(mm,jj)-UCA_nor))*exp(-1i*2*pi*d_ULA_nor*(Nr_ULA_index-1)*cos(DOA_Theta(mm,jj)))*sinc((l-1)-delay(mm,jj))*(-1i*2*pi*R_nor*sin(DOA_Theta(mm,jj))*(-sin(DOA_Phi(mm,jj))));
            end
        end
    end
    Br_angle_Phi_tmp1=cell(1,Nt);
    for jj = 1 : Nt
        Br_angle_Phi_tmp1{1,jj}=Br_angle_Phi_tmp(:,:,jj);
    end
    Br_angle_Phi=blkdiag(Br_angle_Phi_tmp1{:});
end

