clear all;
%close all;
clc; 


%% Declear

Nt      = 1;                    % number of transmit antennas

Nr_UCA  = 16;                   % UCA
Nr_ULA  = 4;                    % UCyA

Nr      = Nr_UCA * Nr_ULA;      % ULA

L       = 4;                    % channel order
N       = 12;
M       = 20;
Pxp     = 1;
K       = 64;                  % OFDM subcarriers
F       = dftmtx(K);
FL      = F(:,1:L);
fc      = 28e9;                 % Carrier Freq (mmWave)
K_r     = db2mag(normrnd(7, 4));
abs_delay = 6 * 10^(-9);
rt      = 3;
gamma_z = db2mag(3);
AOA_LOS = 45;
ZOA_LOS = 45;
c_DS    = 5;
c_ASA   = 17;
c_ZSA   = 7;


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

mag_f = {};
phase_f = {};
delay_f = {};
ZOA_f = {};
AOA_f = {};

CRB_op_ULA_f         = [];
CRB_op_ULA_spec_f    = [];
CRB_op_UCyA_f        = [];
CRB_op_UCyA_spec_f   = [];

for t = 1:10
    tic
    fprintf('Working at: %d iters.\n', t);
    %% Channel generation
    mag         = zeros(N, M, Nt);      % Path loss
    phase       = zeros(N, M, Nt);      % Sub-path phase
    delay       = zeros(N, M, Nt);      % TC + sub-path delay
    ZOA         = zeros(N, M, Nt);      % Sub-path phase offset Theta
    AOA         = zeros(N, M, Nt);      % Sub-path phase offset AOA
    
    for nt = 1 : Nt
        phase_n = unifrnd(-pi, pi, [N, 1]);
        DS      = gen_DS(fc);
        delay_n = gen_cluster_delays(rt, DS, N);
        mag_n   = gen_cluster_powers(delay_n, rt, DS, gamma_z, N);
        ASA     = gen_ASA(fc);
        AOA_n   = gen_cluster_AOAs(ASA, mag_n, K_r, AOA_LOS, N);
        ZSA     = gen_ZSA(fc);
        ZOA_n   = gen_cluster_ZOAs(ZSA, mag_n, K_r, ZOA_LOS, N);
        
        for nn = 1 : N
            for mm = 1 : M
                mag(nn, mm, nt) = mag_n(nn) / M;
                phase(nn, mm, nt) = phase_n(nn);
                switch mm 
                    case {1, 2, 3, 4, 5, 6, 7, 8, 19, 20}
                        delay(nn, mm, nt) = delay_n(nn) + abs_delay;
                    case {9, 10, 11, 12, 17, 18}
                        delay(nn, mm, nt) = delay_n(nn) + 1.28 * c_DS * 10^-9 + abs_delay;
                    case {13, 14, 15, 16}
                        delay(nn, mm, nt) = delay_n(nn) + 2.56 * c_DS * 10^-9 + abs_delay;
                end
                
                switch mm 
                    case {1, 2}
                        offset = 0.0447;
                    case {3, 4}
                        offset = 0.1413;
                    case {5, 6}
                        offset = 0.2492;
                    case {7, 8}
                        offset = 0.3715;
                    case {9, 10}
                        offset = 0.5129;
                    case {11, 12}
                        offset = 0.6797;
                    case {13, 14}
                        offset = 0.8844;
                    case {15, 16}
                        offset = 1.1481;
                    case {17, 18}
                        offset = 1.5195;
                    case {19, 20}
                        offset = 2.1551;
                end
                rand_state_AOA = 2*randi(2) - 3;
                AOA(nn, mm, nt) = (AOA_n(nn) + c_ASA* rand_state_AOA * offset) / 2*pi;
                rand_state_ZOA = 2*randi(2) - 3;
                ZOA(nn, mm, nt) = (ZOA_n(nn) + c_ZSA* rand_state_ZOA * offset) / 2*pi;
            end
        end
    end     
    mag_f{t} = mag;
    phase_f{t} = phase;
    delay_f{t} = delay;
    ZOA_f{t} = ZOA;
    AOA_f{t} = AOA;
    
                
    %% Derivative
    dev_h_mag_ULA       = [];
    dev_h_phase_ULA     = [];
    dev_h_delay_ULA     = [];
    dev_h_ZOA_ULA       = [];
    dev_h_AOA_ULA       = [];
    
    dev_h_mag_UCyA      = [];
    dev_h_phase_UCyA    = [];
    dev_h_delay_UCyA    = [];
    dev_h_ZOA_UCyA      = [];
    dev_h_AOA_UCyA      = [];
    
    for Nr_index=1:Nr
        Br_mag                  = SEMI_spec_chan_derive_mag(        mag, phase, delay, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_mag_ULA           = [dev_h_mag_ULA; transpose(Br_mag)];

        Br_phase                = SEMI_spec_chan_derive_phase(      mag, phase, delay, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_phase_ULA         = [dev_h_phase_ULA; transpose(Br_phase)];

        Br_delay                = SEMI_spec_chan_derive_delay(      mag, phase, delay, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_delay_ULA         = [dev_h_delay_ULA; transpose(Br_delay)];
        
        Br_angle_ZOA            = SEMI_spec_chan_derive_ZOA(        mag, phase, delay, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_ZOA_ULA           = [dev_h_ZOA_ULA; transpose(Br_angle_ZOA)];

        Br_angle_AOA            = SEMI_spec_chan_derive_AOA(        mag, phase, delay, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), L, N, M, Nt);
        dev_h_AOA_ULA           = [dev_h_AOA_ULA; transpose(Br_angle_AOA)];
    end
    
    for Nr_ULA_index=1:Nr_ULA
        for Nr_UCA_index=1:Nr_UCA
            Br_mag              = SEMI_spec_chan_derive_mag(        mag, phase, delay, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_mag_UCyA      = [dev_h_mag_UCyA; transpose(Br_mag)];

            Br_phase            = SEMI_spec_chan_derive_phase(      mag, phase, delay, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_phase_UCyA    = [dev_h_phase_UCyA; transpose(Br_phase)];

            Br_delay            = SEMI_spec_chan_derive_delay(      mag, phase, delay, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_delay_UCyA    = [dev_h_delay_UCyA; transpose(Br_delay)];
            
            Br_angle_ZOA        = SEMI_spec_chan_derive_ZOA(        mag, phase, delay, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_ZOA_UCyA      = [dev_h_ZOA_UCyA; transpose(Br_angle_ZOA)];

            Br_angle_AOA        = SEMI_spec_chan_derive_AOA(        mag, phase, delay, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), L, N, M, Nt);
            dev_h_AOA_UCyA      = [dev_h_AOA_UCyA; transpose(Br_angle_AOA)];
        end
    end

    %% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters
    G_ULA   = [dev_h_mag_ULA,  dev_h_phase_ULA,  dev_h_delay_ULA,  dev_h_ZOA_ULA,  dev_h_AOA_ULA]; 
    G_UCyA  = [dev_h_mag_UCyA, dev_h_phase_UCyA, dev_h_delay_UCyA, dev_h_ZOA_UCyA, dev_h_AOA_UCyA]; 
    
    %% CRB
    N_total = K;
    N_pilot = N_total/4;
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
        try
            CRB_op_ULA_spec(snr_i)  = abs(trace(pinv(Iop_spec_ULA)));
        catch ME
            t = t - 1;
            continue
        end
        
        Iop_spec_UCyA           = G_UCyA*G_UCyA'*Iop_full*G_UCyA*G_UCyA';
        try
            CRB_op_UCyA_spec(snr_i) = abs(trace(pinv(Iop_spec_UCyA)));
        catch ME
            t = t - 1;
            continue
        end
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