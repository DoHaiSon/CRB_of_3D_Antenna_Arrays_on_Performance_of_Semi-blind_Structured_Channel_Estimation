function Br_angle_AOA= SEMI_spec_chan_derive_AOA(gamma, AOA, ZOA, position_elements_nor, d_H, Nt)
%Nt   % number of transmit antennas
%Nr   % number of receive antennas
%L    % channel order
%N    % number of TCs
%M    % number of sub-paths per TC

    Br_angle_AOA_tmp = zeros(d_H,Nt);
    for jj = 1 : Nt
        for nn = 1 : d_H
            r_x = sin(ZOA(nn,jj)) * cos(AOA(nn,jj));
            r_y = sin(ZOA(nn,jj)) * sin(AOA(nn,jj));
            r_z = cos(ZOA(nn,jj));
            Br_angle_AOA_tmp(nn,jj)=gamma(nn,jj) ...
                * exp(-1i*2*pi*(position_elements_nor(1)*r_x + position_elements_nor(2)*r_y + position_elements_nor(3)*r_z)) ...
                * (-1i*2*pi)*(-position_elements_nor(1) * sin(ZOA(nn,jj)) * sin(AOA(nn,jj)) + ...
                + position_elements_nor(2) * sin(ZOA(nn,jj)) * cos(AOA(nn,jj))  + position_elements_nor(3) * r_z);
        end
    end
    
    Br_angle_AOA_tmp1 = cell(1,Nt);
    for nt = 1 : Nt
        Br_angle_AOA_tmp1{1, nt} = blkdiag(Br_angle_AOA_tmp(:, nt));
    end
    
    Br_angle_AOA = blkdiag(Br_angle_AOA_tmp1{:});
end

