clear all;
%close all;
clc; 


%% Declear

Nt      = 2;                    % number of transmit antennas

Nr_UCA  = 16;                   % UCA
Nr_ULA  = 4;                    % UCyA

Nr      = Nr_UCA * Nr_ULA;      % ULA

L       = 4;                    % channel order
% N       = ceil(unifrnd(1, 10)); % Number of TCs
% M       = ceil(unifrnd(1, 10)); % Number of sub-path per TC  
N       = 4;
M       = 4;
Pxp     = 1;
K       = 256;                  % OFDM subcarriers
F       = dftmtx(K);
FL      = F(:,1:L);
sigmax2 = 1;
fc      = 28e9;                 % Carrier Freq (mmWave)

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

CRB_op_ULA_f         = [];
CRB_op_ULA_spec_f    = [];
CRB_op_UCyA_f        = [];
CRB_op_UCyA_spec_f   = [];

for t = 1:2
    tic
    fprintf('Working at: %d iters.\n', t);
    %% Channel generation
    mag         = zeros(N, M, Nt);      % Path loss
    phase       = zeros(N, M, Nt);      % Sub-path phase
    delay       = zeros(N, M, Nt);      % TC + sub-path delay
    DOA_Theta   = zeros(N, M, Nt);      % Sub-path phase offset Theta
    DOA_Phi     = zeros(N, M, Nt);      % Sub-path phase offset Phi
    
    for nt = 1 : Nt
        phase_nt =  unifrnd(0, 2*pi, [N, 1]);
%         delay_nt = 
        for nn = 1 : N
            for mm = 1 : M
                phase(nt, nn, mm) = phase_nt(nn);
                               
            end
        end
    end
    mag         = gen_pathloss(fc, 100, 1.2, 1.8, N, M, Nt);
    delay       = exprnd(1/0.043, [N, M, Nt])  * 10^(-9);
    DOA_Theta   = normrnd(0, 19.3, [N, M, Nt]) * (pi/180);
    DOA_Phi     = normrnd(0, 11.3, [N, M, Nt]) * (pi/180);       
                
    %% Derivative
    dev_h_mag_ULA           = [];
    dev_h_phase_ULA         = [];
    dev_h_delay_ULA         = [];
    dev_h_angle_Theta_ULA   = [];
    dev_h_angle_Phi_ULA     = [];
    
    dev_h_mag_UCyA          = [];
    dev_h_phase_UCyA        = [];
    dev_h_delay_UCyA        = [];
    dev_h_angle_Theta_UCyA  = [];
    dev_h_angle_Phi_UCyA    = [];
    
    for Nr_index=1:Nr
        Br_mag                  = SEMI_spec_chan_derive_mag(        mag, phase, delay, DOA_Phi, DOA_Theta, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_mag_ULA           = [dev_h_mag_ULA; transpose(Br_mag)];

        Br_phase                = SEMI_spec_chan_derive_phase(      mag, phase, delay, DOA_Phi, DOA_Theta, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_phase_ULA         = [dev_h_phase_ULA; transpose(Br_phase)];

        Br_delay                = SEMI_spec_chan_derive_delay(      mag, phase, delay, DOA_Phi, DOA_Theta, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_delay_ULA         = [dev_h_delay_ULA; transpose(Br_delay)];
        
        Br_angle_Theta          = SEMI_spec_chan_derive_angle_Theta(mag, phase, delay, DOA_Phi, DOA_Theta, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_angle_Theta_ULA   = [dev_h_angle_Theta_ULA; transpose(Br_angle_Theta)];

        Br_angle_Phi            = SEMI_spec_chan_derive_angle_Phi(  mag, phase, delay, DOA_Phi, DOA_Theta, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_angle_Phi_ULA     = [dev_h_angle_Phi_ULA; transpose(Br_angle_Phi)];
    end
    
    for Nr_ULA_index=1:Nr_ULA
        for Nr_UCA_index=1:Nr_UCA
            Br_mag              = SEMI_spec_chan_derive_mag(        mag, phase, delay, DOA_Phi, DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_mag_UCyA      = [dev_h_mag_UCyA; transpose(Br_mag)];

            Br_phase            = SEMI_spec_chan_derive_phase(      mag, phase, delay, DOA_Phi, DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_phase_UCyA    = [dev_h_phase_UCyA; transpose(Br_phase)];

            Br_delay            = SEMI_spec_chan_derive_delay(      mag, phase, delay, DOA_Phi, DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_delay_UCyA    = [dev_h_delay_UCyA; transpose(Br_delay)];
            
            Br_angle_Theta      = SEMI_spec_chan_derive_angle_Theta(mag, phase, delay, DOA_Phi, DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_angle_Theta_UCyA = [dev_h_angle_Theta_UCyA; transpose(Br_angle_Theta)];

            Br_angle_Phi        = SEMI_spec_chan_derive_angle_Phi(  mag, phase, delay, DOA_Phi, DOA_Theta, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_angle_Phi_UCyA   = [dev_h_angle_Phi_UCyA; transpose(Br_angle_Phi)];
        end
    end

    %% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters
    G_ULA   = [dev_h_mag_ULA,  dev_h_phase_ULA,  dev_h_delay_ULA,  dev_h_angle_Theta_ULA,  dev_h_angle_Phi_ULA]; 
    G_UCyA  = [dev_h_mag_UCyA, dev_h_phase_UCyA, dev_h_delay_UCyA, dev_h_angle_Theta_UCyA, dev_h_angle_Phi_UCyA]; 
    
    %% CRB
    N_total = 64;
    N_pilot = 64/4;
    N_data  = N_total-N_pilot;
    %============================================
    SNR = -10:5:20;
    for snr_i = 1 : length(SNR)
        sigmav2 = 10^(-SNR(snr_i)/10);
    %============================================
        %Only Pilot    
        X_nga = kron(eye(Nr * N),X);
    %============================================
        Iop                     = X_nga' * X_nga / sigmav2;
        Iop_full                = N_pilot * Iop;
        CRB_op(snr_i)           = abs(trace(pinv(Iop_full)));
    %============================================
        %Only Pilot Specular   
        Iop_spec_ULA            = G_ULA*G_ULA'*Iop_full*G_ULA*G_ULA';
        CRB_op_ULA_spec(snr_i)  = abs(trace(pinv(Iop_spec_ULA)));
        
        Iop_spec_UCyA           = G_UCyA*G_UCyA'*Iop_full*G_UCyA*G_UCyA';
        CRB_op_UCyA_spec(snr_i) = abs(trace(pinv(Iop_spec_UCyA)));
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
figure()
semilogy(SNR, CRB_op_ULA, '-b');
hold on;
semilogy(SNR, CRB_op_UCyA, '*');
semilogy(SNR, CRB_op_ULA_spec, '-g o');
semilogy(SNR, CRB_op_UCyA_spec, '-r d');

grid on;
ylabel('Average CRB (dB)');
xlabel('SNR (dB)');
legend('ULA: unst_{P}', 'UCyA: unst_{P}', 'ULA: stru_{P}', 'UCyA: stru_{P}');
hold on;