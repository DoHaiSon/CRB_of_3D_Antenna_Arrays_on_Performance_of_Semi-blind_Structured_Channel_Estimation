clear all;
%close all;
clc; 

CRB_op_ULA_f         = [];
CRB_op_ULA_spec_f    = [];
CRB_op_UCyA_f        = [];
CRB_op_UCyA_spec_f   = [];

for t = 1:100
    tic
    fprintf('Working at: %d iters.\n', t);
    
    %% Declear

    Nt      = 2;         % number of transmit antennas

    Nr_UCA  = 16;        % UCA
    Nr_ULA  = 4;         % UCyA

    Nr      = Nr_UCA * Nr_ULA;        % ULA

    L       = 4;         % channel order
    for M = 1:5          % Number of multipaths  
        Pxp     = 1;
        K       = 64;        % OFDM subcarriers
        F       = dftmtx(K);
        FL      = F(:,1:L);
        sigmax2 = 1;

        %% Generate position of elements in arrays
        d_ULA_nor   = 0.5;
        d_UCA_nor   = 0.5;
        R_nor       = 0.5 * d_UCA_nor/sin(pi/Nr_UCA);

        ULA_elements_nor    = zeros(3, 1, Nr);
        UCyA_elements_nor   = zeros(3, Nr_ULA, Nr_UCA);

        for Nr_index=1:Nr
            ULA_elements_nor(1, 1, Nr_index) = (Nr_index-1) * d_UCA_nor;
            ULA_elements_nor(2, 1, Nr_index) = 0;
            ULA_elements_nor(3, 1, Nr_index) = 0;
        end

        for Nr_ULA_index=1:Nr_ULA
            for Nr_UCA_index=1:Nr_UCA
                UCyA_elements_nor(1, Nr_ULA_index, Nr_UCA_index) = R_nor * sin((Nr_UCA_index-1)*(2*pi/Nr_UCA)) ;         % x
                UCyA_elements_nor(2, Nr_ULA_index, Nr_UCA_index) = R_nor * cos((Nr_UCA_index-1)*(2*pi/Nr_UCA)) ;         % y
                UCyA_elements_nor(3, Nr_ULA_index, Nr_UCA_index) = (Nr_ULA_index-1) * d_ULA_nor;                         % z
            end
        end

        %% Signal Generation
        % we use the Zadoff-Chu sequences
        U = 1:2:7;
        ZC_p = [];
        for u = 1 : Nt
            for k = 1 : K
                ZC(k,u) = sqrt(Pxp) * exp( ( -1i * pi * U(u) * (k-1)^2 ) / K );
            end
            ZC_p = [ZC_p; ZC(:,u)];
        end

        X = [];
        for ii = 1 : Nt
            X  = [X diag(ZC(:,ii))*FL];
        end

        %% Channel generation 
        % Fading, delay, DOA matrix of size(M,Nt), M - the number of multipath
        % Assume that real and imaginary parts are independent
        % 1/sqrt(2) makes the variaonce to be 1
        fading      = 1/sqrt(2)*(rand(M, Nt)  + 1i*rand(M, Nt));
        delay       = exprnd(1/0.043, M, Nt)  * 10^(-9);
        DOA_Phi     = normrnd(0, 11.3, M, Nt) * (pi/180);
        DOA_Theta   = normrnd(0, 19.3, M, Nt) * (pi/180);

        %% Derivative
        dev_h_fading_ULA        = [];
        dev_h_conj_fading_ULA   = [];
        dev_h_delay_ULA         = [];
        dev_h_angle_Phi_ULA     = [];
        dev_h_angle_Theta_ULA   = [];

        dev_h_fading_UCyA       = [];
        dev_h_conj_fading_UCyA  = [];
        dev_h_delay_UCyA        = [];
        dev_h_angle_Phi_UCyA    = [];
        dev_h_angle_Theta_UCyA  = [];

        for Nr_index=1:Nr
            Br_fading               = SEMI_spec_chan_derive_fading(fading,delay,DOA_Phi,DOA_Theta, ULA_elements_nor(:, 1, Nr_index),L,M,Nt);
            dev_h_fading_ULA        = [dev_h_fading_ULA; transpose(Br_fading)];

            Br_conj_fading          = SEMI_spec_chan_derive_conj_fading(fading,delay,DOA_Phi,DOA_Theta, ULA_elements_nor(:, 1, Nr_index),L,M,Nt);
            dev_h_conj_fading_ULA   = [dev_h_conj_fading_ULA; transpose(Br_conj_fading)];

            Br_delay                = SEMI_spec_chan_derive_delay(fading,delay,DOA_Phi,DOA_Theta,ULA_elements_nor(:, 1, Nr_index),L,M,Nt);
            dev_h_delay_ULA         = [dev_h_delay_ULA; transpose(Br_delay)];

            Br_angle_Phi            = SEMI_spec_chan_derive_angle_Phi(fading,delay,DOA_Phi,DOA_Theta,ULA_elements_nor(:, 1, Nr_index),L,M,Nt);
            dev_h_angle_Phi_ULA     = [dev_h_angle_Phi_ULA; transpose(Br_angle_Phi)];

            Br_angle_Theta          = SEMI_spec_chan_derive_angle_Theta(fading,delay,DOA_Phi,DOA_Theta,ULA_elements_nor(:, 1, Nr_index),L,M,Nt);
            dev_h_angle_Theta_ULA   = [dev_h_angle_Theta_ULA; transpose(Br_angle_Theta)];
        end

        for Nr_ULA_index=1:Nr_ULA
            for Nr_UCA_index=1:Nr_UCA
                Br_fading           = SEMI_spec_chan_derive_fading(fading,delay,DOA_Phi,DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
                dev_h_fading_UCyA   = [dev_h_fading_UCyA; transpose(Br_fading)];

                Br_conj_fading      = SEMI_spec_chan_derive_conj_fading(fading,delay,DOA_Phi,DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
                dev_h_conj_fading_UCyA = [dev_h_conj_fading_UCyA; transpose(Br_conj_fading)];

                Br_delay            = SEMI_spec_chan_derive_delay(fading,delay,DOA_Phi,DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
                dev_h_delay_UCyA    = [dev_h_delay_UCyA; transpose(Br_delay)];

                Br_angle_Phi        = SEMI_spec_chan_derive_angle_Phi(fading,delay,DOA_Phi,DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
                dev_h_angle_Phi_UCyA   = [dev_h_angle_Phi_UCyA; transpose(Br_angle_Phi)];

                Br_angle_Theta      = SEMI_spec_chan_derive_angle_Theta(fading,delay,DOA_Phi,DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
                dev_h_angle_Theta_UCyA = [dev_h_angle_Theta_UCyA; transpose(Br_angle_Theta)];
            end
        end

        %% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters
        G_ULA   = [dev_h_fading_ULA,  dev_h_conj_fading_ULA,  dev_h_delay_ULA,  dev_h_angle_Theta_ULA,  dev_h_angle_Phi_ULA]; 
        G_UCyA  = [dev_h_fading_UCyA, dev_h_conj_fading_UCyA, dev_h_delay_UCyA, dev_h_angle_Theta_UCyA, dev_h_angle_Phi_UCyA]; 

        %% CRB
        N_total = 64;
        N_pilot = 64/4;
        N_data  = N_total-N_pilot;
        %============================================
        SNR = 5;
        sigmav2 = 10^(-SNR/10);
    %============================================
        %Only Pilot    
        X_nga = kron(eye(Nr),X);
    %============================================
        Iop                     = X_nga' * X_nga / sigmav2;
        Iop_full                = N_pilot * Iop;
        CRB_op(M)               = abs(trace(pinv(Iop_full)));
    %============================================
        %Only Pilot Specular   
        Iop_spec_ULA            = G_ULA*G_ULA'*Iop_full*G_ULA*G_ULA';
        CRB_op_ULA_spec(M)      = abs(trace(pinv(Iop_spec_ULA)));

        Iop_spec_UCyA           = G_UCyA*G_UCyA'*Iop_full*G_UCyA*G_UCyA';
        CRB_op_UCyA_spec(M)     = abs(trace(pinv(Iop_spec_UCyA)));
    end
    
    CRB_op_ULA_f        = [CRB_op_ULA_f; CRB_op];
    CRB_op_UCyA_f       = [CRB_op_UCyA_f; CRB_op];
    CRB_op_ULA_spec_f   = [CRB_op_ULA_spec_f; CRB_op_ULA_spec];
    CRB_op_UCyA_spec_f  = [CRB_op_UCyA_spec_f; CRB_op_UCyA_spec];
    toc
end

CRB_op_ULA      = mean(CRB_op_ULA_f);
CRB_op_UCyA     = mean(CRB_op_UCyA_f);
CRB_op_ULA_spec = mean(CRB_op_ULA_spec_f);
CRB_op_UCyA_spec= mean(CRB_op_UCyA_spec_f);

%% figure
M = 1:5;
figure()
semilogy(M, CRB_op_ULA, '-k');
hold on;
semilogy(M, CRB_op_UCyA, '*');
semilogy(M, CRB_op_ULA_spec, '-g o');
semilogy(M, CRB_op_UCyA_spec, '-r d');

grid on;
ylabel('Average CRB (dB)');
xlabel('M');
legend('ULA: unst', 'UCyA: unst', 'ULA: stru', 'UCyA: stru');
hold on;