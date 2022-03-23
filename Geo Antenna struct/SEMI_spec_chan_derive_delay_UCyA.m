function Br_delay= SEMI_spec_chan_derive_delay_UCyA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor,L,M,Nt)

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
Br_delay_tmp = zeros(M,L,Nt);
for jj = 1 : Nt
    for mm = 1 : M
        for l = 1 : L
            r_x = sin(DOA_Theta(mm,jj)) * cos(DOA_Phi(mm,jj));
            r_y = sin(DOA_Theta(mm,jj)) * sin(DOA_Phi(mm,jj));
            r_z = cos(DOA_Theta(mm,jj));
            Br_delay_tmp(mm,l,jj)=fading(mm,jj)*((sin((l-1)-delay(mm,jj))/(((l-1)-delay(mm,jj))^2))-(cos((l-1)-delay(mm,jj))/((l-1)-delay(mm,jj)))) ...
                * exp(-1i*2*pi*(position_elements_nor(1)*r_x + position_elements_nor(2)*r_y + position_elements_nor(3)*r_z));
        end
    end
end
Br_delay_tmp1=cell(1,Nt);
for jj = 1 : Nt
Br_delay_tmp1{1,jj}=Br_delay_tmp(:,:,jj);
end
Br_delay=blkdiag(Br_delay_tmp1{:});
end

