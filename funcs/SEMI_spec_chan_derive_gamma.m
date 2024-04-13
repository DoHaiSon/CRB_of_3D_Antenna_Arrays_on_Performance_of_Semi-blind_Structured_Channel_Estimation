function Br_mag= SEMI_spec_chan_derive_gamma(gamma, AOA, ZOA, position_elements_nor, d_H, Nt)
%Nt   % number of transmit antennas
%Nr   % number of receive antennas
%d_H  % number of channels

    Br_mag_tmp = zeros(d_H,Nt);
    for jj = 1 : Nt
        for nn = 1: d_H
            r_x = sin(ZOA(nn,jj)) * cos(AOA(nn,jj));
            r_y = sin(ZOA(nn,jj)) * sin(AOA(nn,jj));
            r_z = cos(ZOA(nn,jj));
            Br_mag_tmp(nn,jj)= 1/2 * (1-1i) ...
                * exp(-1i*2*pi*(position_elements_nor(1)*r_x + position_elements_nor(2)*r_y + position_elements_nor(3)*r_z));
        end
    end
    
    Br_mag_tmp1 = cell(1,Nt);
    for nt = 1 : Nt
        Br_mag_tmp1{1, nt} = blkdiag(Br_mag_tmp(:,nt));
    end
    
    Br_mag = blkdiag(Br_mag_tmp1{:});
end