clear all;
%close all;
clc; 


%% Declear

Nt      = 2;         % number of transmit antennas
Nr_UCA  = 8;         % number of receive antennas of UCA
Nr_ULA  = 2;         % number of receive antennas of ULA
L       = 4;         % channel order
M       = 2;         % Number of multipaths   
Pxp     = 0.3;
K       = 64;        % OFDM subcarriers
F       = dftmtx(K);
FL      = F(:,1:L);
sigmax2 = 1;
fc      = 0;

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

fading      = rand(M,Nt);
delay       = rand(M,Nt);
DOA_Phi     = rand(M,Nt) * pi;
DOA_Theta   = rand(M,Nt) * pi;
d_ULA_nor   = 0.5;
d_UCA_nor   = 0.5;
R_nor       = 0.5 * d_UCA_nor/sin(pi/Nr_UCA);

ULA_nor                 = zeros(Nr_ULA, Nr_UCA);
position_elements_nor   = zeros(3, Nr_ULA, Nr_UCA);

for Nr_ULA_index=1:Nr_ULA
    for Nr_UCA_index=1:Nr_UCA
        position_elements_nor(1, Nr_ULA_index, Nr_UCA_index) = (Nr_UCA_index-1) * d_UCA_nor;
        position_elements_nor(2, Nr_ULA_index, Nr_UCA_index) = 0;
        position_elements_nor(3, Nr_ULA_index, Nr_UCA_index) = 0;
%         position_elements_nor(1, Nr_ULA_index, Nr_UCA_index) = R_nor * sin((Nr_UCA_index-1)*(2*pi/Nr_UCA)) ;         % x
%         position_elements_nor(2, Nr_ULA_index, Nr_UCA_index) = R_nor * cos((Nr_UCA_index-1)*(2*pi/Nr_UCA)) ;         % y
%         position_elements_nor(3, Nr_ULA_index, Nr_UCA_index) = (Nr_ULA_index-1) * d_ULA_nor;                         % z
%         ULA_nor(Nr_ULA_index, Nr_UCA_index) = (Nr_UCA_index-1)*2*pi/Nr_UCA + rot_nor * (Nr_ULA_index-1)*2*pi/Nr_UCA;    % polar axis
    end
end

%% Derivative

dev_h_fading        = [];
dev_h_delay         = [];
dev_h_angle_Phi     = [];
dev_h_angle_Theta   = [];

for Nr_ULA_index=1:Nr_ULA
    for Nr_UCA_index=1:Nr_UCA
        Br_fading         = SEMI_spec_chan_derive_fading_UCyA(fading,delay,DOA_Phi,DOA_Theta, position_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt,fc);
        dev_h_fading      = [dev_h_fading; transpose(Br_fading)];

        Br_delay          = SEMI_spec_chan_derive_delay_UCyA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt,fc);
        dev_h_delay       = [dev_h_delay; transpose(Br_delay)];

        Br_angle_Phi      = SEMI_spec_chan_derive_angle_Phi_UCyA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt,fc);
        dev_h_angle_Phi   = [dev_h_angle_Phi; transpose(Br_angle_Phi)];

        Br_angle_Theta    = SEMI_spec_chan_derive_angle_Theta_UCyA(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor(:, Nr_ULA_index, Nr_UCA_index),L,M,Nt,fc);
        dev_h_angle_Theta = [dev_h_angle_Theta; transpose(Br_angle_Theta)];
    end
end

%% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters
G = [dev_h_fading, dev_h_delay, dev_h_angle_Theta, dev_h_angle_Phi]; 

%% ------------------------------------------------------------------------
 
X = [];
for ii = 1 : Nt
    X  = [X diag(ZC(:,ii))*FL];
end

[H, h_true] = gen_chan_specular(fading,delay,DOA_Phi,DOA_Theta,position_elements_nor,Nr_UCA,Nr_ULA,L,Nt,fc);


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

%figure
figure()
semilogy(SNR,CRB_op,'-b +');
hold on;
semilogy(SNR,CRB_op_spec,'-r *');

semilogy(SNR,CRB_SB,'-g +');
hold on;
semilogy(SNR,CRB_SB_spec,'-k *');

grid on;
ylabel('Normalized CRB');
xlabel('SNR(dB)');
legend('normal_{OP}','spec_{OP}','normal_{SB}','spec_{SB}');
title(' ');
hold on;