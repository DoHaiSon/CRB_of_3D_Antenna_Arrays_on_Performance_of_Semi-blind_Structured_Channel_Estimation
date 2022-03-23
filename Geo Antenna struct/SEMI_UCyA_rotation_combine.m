clear all;
%close all;
clc; 


%% Declear

Nt      = 2;         % number of transmit antennas
Nr_UCA  = 4;         % number of receive antennas of UCA
Nr_ULA  = 2;         % number of receive antennas of ULA
L       = 4;         % channel order
M       = 2;         % Number of multipaths (assumption: M  = L)   
Pxp     = 0.3;
K       = 64;        % OFDM subcarriers
F       = dftmtx(K);
FL      = F(:,1:L);
sigmax2 = 1;


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

%% Channel generation 
% Fading, delay, DOA matrix of size(M,Nt), M - the number of multipath

fading      = [0.8, 0.6, 0.4, 0.2; 0.9, 0.7, 0.5, 0.3];
delay       = [1, 2, 3, 4; 2, 3, 4, 5]*10^(-9);
DOA_Phi     = [0.3, 0.4, 0.25, 0.6; 0.7, 0.85, 0.43, 0.66] * pi;
DOA_Theta   = [pi/2, pi/4, pi/6, pi/8; pi/3, pi/5, pi/7, pi/9];
d_ULA_nor   = 0.5;
d_UCA_nor   = 0.5;
R_nor       = 0.5 * d_UCA_nor/sin(pi/Nr_UCA);

step = 1:10;

for rotation = 1:length(step)
    rot_nor = (rotation - 1) / length(step);
    fprintf('\n--------------------------------------------------\nWorking at rotation: %f.\n\n', rot_nor);
    
    ULA_nor = zeros(Nr_ULA, Nr_UCA);
    position_elements_nor = zeros(3, Nr_ULA, Nr_UCA);
    for Nr_ULA_index=1:Nr_ULA
        for Nr_UCA_index=1:Nr_UCA
            position_elements_nor(1, Nr_ULA_index, Nr_UCA_index) = R_nor * sin((Nr_UCA_index-1)*(2*pi/Nr_UCA) + rot_nor * (2*pi/Nr_UCA)) ;         % x
            position_elements_nor(2, Nr_ULA_index, Nr_UCA_index) = R_nor * cos((Nr_UCA_index-1)*(2*pi/Nr_UCA) + rot_nor * (2*pi/Nr_UCA)) ;         % y
            position_elements_nor(3, Nr_ULA_index, Nr_UCA_index) = (Nr_ULA_index-1) * d_ULA_nor;                                                   % z
        end
    end

    %% Derivative

    dev_h_fading        = [];
    dev_h_delay         = [];
    dev_h_angle_Phi     = [];
    dev_h_angle_Theta   = [];

    for Nr_ULA_index=1:Nr_ULA
        for Nr_UCA_index=1:Nr_UCA
            Br_fading         = SEMI_spec_chan_derive_fading_UCyA(fading,delay,DOA_Phi,DOA_Theta, position_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
            dev_h_fading      = [dev_h_fading; transpose(Br_fading)];
            
            Br_delay          = SEMI_spec_chan_derive_delay_UCyA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
            dev_h_delay       = [dev_h_delay; transpose(Br_delay)];
            
            Br_angle_Theta    = SEMI_spec_chan_derive_angle_Theta_UCyA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
            dev_h_angle_Theta = [dev_h_angle_Theta; transpose(Br_angle_Theta)];
            
            Br_angle_Phi      = SEMI_spec_chan_derive_angle_Phi_UCyA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt);
            dev_h_angle_Phi   = [dev_h_angle_Phi; transpose(Br_angle_Phi)];
        end
    end

    %% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters
    G = [dev_h_fading, dev_h_delay, dev_h_angle_Theta, dev_h_angle_Phi]; 
    %% ------------------------------------------------------------------------

    X = [];
    for ii = 1 : Nt
        X  = [X diag(ZC(:,ii))*FL];
    end

    [H, h_true] = gen_chan_specular_rotation(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor,Nr_UCA,Nr_ULA,L,Nt);


    %% LAMBDA
    LAMBDA  = [];
    for jj = 1 : Nt
        lambda_j =[];
        for r = 1 : Nr_ULA*Nr_UCA
            h_rj       = transpose(H(r,:,jj));
            lambda_rj  = diag(FL*h_rj);
            lambda_j   = [lambda_j; lambda_rj];
        end
        LAMBDA = [LAMBDA lambda_j];
    end


    %% Partial derivative of LAMBDA w.r.t. h_i
    partial_LAMBDA  = cell(1, L);
    for ll = 1 : L
        partial_LAMBDA_ll = [];
        for jj = 1 : Nt
            lambda_jj =[];
            for r = 1 : Nr_ULA*Nr_UCA
                lambda_rj_ll = diag(FL(:,ll));
                lambda_jj    = [lambda_jj; lambda_rj_ll];
            end
            partial_LAMBDA_ll = [partial_LAMBDA_ll lambda_jj];
        end
        partial_LAMBDA{1,ll} = partial_LAMBDA_ll;
    end

    N_total = 4;
    N_pilot = 2;
    N_data  = N_total-N_pilot;
    %============================================
    SNR = -10:5:30;
    for snr_i = 1 : length(SNR)
        tic
        fprintf('Working at SNR: %d dB.\n', SNR(snr_i));

        sigmav2 = 10^(-SNR(snr_i)/10);
    %============================================
        %Only Pilot    
        X_nga=kron(eye(Nr_ULA*Nr_UCA),X);
    %============================================
        Iop           = X_nga'*X_nga / sigmav2;
        Iop_full      = N_total * Iop;
        CRB_op(snr_i) = abs(trace(pinv(Iop_full)));
    %============================================
        %Only Pilot Specular   
        Iop_spec = G*G'*Iop_full*G*G';
        CRB_op_spec(snr_i) = abs(trace(pinv(Iop_spec)));
    %============================================
    %% SemiBlind
        Cyy      = sigmax2 * LAMBDA * LAMBDA'  + sigmav2 * eye(K*Nr_ULA*Nr_UCA);
        Cyy_inv  = pinv(Cyy);   

        for ii = 1 : L
            partial_Cyy_hii = sigmax2 * LAMBDA * partial_LAMBDA{1,ii}';
            for jj = ii : L
                partial_Cyy_hjj = sigmax2 * LAMBDA * partial_LAMBDA{1,jj}';
                % Slepian-Bangs Formula
                I_D(ii,jj) = trace(Cyy_inv * partial_Cyy_hii * Cyy_inv * partial_Cyy_hjj);
                I_D(ii,jj) = I_D(ii,jj)';
            end
        end
        I_D = triu(repmat(I_D, Nr_ULA*Nr_UCA*Nt, Nr_ULA*Nr_UCA*Nt));
    %============================================
        %Semiblind Normal
        I_SB               = N_data*I_D+N_pilot*Iop;
        CRB_SB_i           = pinv(I_SB);
        CRB_SB(snr_i)      = abs(trace(CRB_SB_i));

        clear I_D;    
    %============================================
        %Semiblind Specular
        I_SB_spec               = G*G'*I_SB*G*G';
        CRB_SB_spec_i           = pinv(I_SB_spec);
        CRB_SB_spec(snr_i)      = abs(trace(CRB_SB_spec_i));   
        toc
    end
    CRB_op_f(rotation, :)       = CRB_op;
    CRB_op_spec_f(rotation, :)  = CRB_op_spec;
    CRB_SB_f(rotation, :)       = CRB_SB;
    CRB_SB_spec_f(rotation, :)  = CRB_SB_spec;
    
end

SNR_f = repmat(SNR, [length(step), 1]);
rot_nor = repmat((1:length(step))', 1, length(SNR));

xlabels = {};
for i = 1:length(step)
    xlabels{end + 1} = mat2str((i-1) / 10);
end
xtick = linspace(0,length(step)-1,length(step));

%figure
figure(1)
grid on;
surf(rot_nor, SNR_f, CRB_op_f);
xlabel('Normalize rotated');
set(gca, 'XTick', xtick);
set(gca, 'XTickLabel', xlabels);
ylabel('SNR(dB)');
zlabel('CRB(dB)')
set(gca,'ZScale','log');
title('Normal OP');

figure(2)
grid on;
surf(rot_nor, SNR_f, CRB_op_spec_f);
xlabel('Normalize rotated');
set(gca, 'XTick', xtick);
set(gca, 'XTickLabel', xlabels);
ylabel('SNR(dB)');
zlabel('CRB(dB)');
set(gca,'ZScale','log');
title('Normal spec');

figure(3)
grid on;
surf(rot_nor, SNR_f, CRB_SB_f);
xlabel('Normalize rotated');
set(gca, 'XTick', xtick);
set(gca, 'XTickLabel', xlabels);
ylabel('SNR(dB)');
zlabel('CRB(dB)');
set(gca,'ZScale','log');
title('Semi-blind OP');

figure(4)
grid on;
surf(rot_nor, SNR_f, CRB_SB_spec_f);
xlabel('Normalize rotated');
set(gca, 'XTick', xtick);
set(gca, 'XTickLabel', xlabels);
ylabel('SNR(dB)');
zlabel('CRB(dB)');
set(gca,'ZScale','log');
title('Semi-blind spec');