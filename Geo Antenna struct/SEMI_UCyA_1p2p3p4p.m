clear all;
%close all;
clc; 


%% Declear

Nt      = 2;         % number of transmit antennas
Nr_UCA  = 2;         % number of receive antennas of UCA
Nr_ULA  = 8;         % number of receive antennas of ULA
L       = 8;         % channel order
M       = 2;         % Number of multipaths (assumption: M  = L)   
Pxp     = 0.3;
K       = 64;        % OFDM subcarriers
F       = dftmtx(K);
FL      = F(:,1:L);


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
%fading = rand(M,Nt)+1i*rand(M,Nt);
%fading = rand(M,Nt);
%fading = 0.1+(0.6-0.1)*rand(M,Nt);
fading=[0.8,0.6,0.4,0.2;0.9,0.7,0.5,0.3];
%delay  = 0.1+(0.15-0.1)*rand(M,Nt)*0.01;
delay=[0.1,0.2,0.3,0.4;0.2,0.3,0.4,0.5]*0.001;
%DOA    = (-0.5+(0.5-(-0.5))*rand(M,Nt))*pi;
DOA_Phi=[pi/2,pi/4,pi/6,pi/8;pi/3,pi/5,pi/7,pi/9];
DOA_Theta=[0.3,0.4,0.25,0.6;0.7,0.85,0.43,0.66]*pi;
%DOA_Theta = (0.1+(1-0.1)*rand(M,Nt))*pi;
d_ULA_nor=0.5;
d_UCA_nor=0.5;
R_nor=0.5*d_UCA_nor/sin(pi/Nr_UCA);
sigmax2=1;

%% Derivative

dev_h_fading_tmp=[];
dev_h_delay_tmp=[];
dev_h_angle_tmp=[];


dev_h_fading=[];
dev_h_delay=[];
dev_h_angle_Phi=[];
dev_h_angle_Theta=[];

for Nr_ULA_index=1:Nr_ULA
    for Nr_UCA_index=1:Nr_UCA
        Br_fading = SEMI_spec_chan_derive_fading_UCyA(fading,delay,DOA_Phi,DOA_Theta,R_nor,d_ULA_nor,Nr_UCA_index,Nr_ULA_index,Nr_UCA,Nr_ULA,L,M,Nt);
        dev_h_fading=[dev_h_fading; transpose(Br_fading)];

        Br_delay = SEMI_spec_chan_derive_delay_UCyA(fading,delay,DOA_Phi,DOA_Theta,R_nor,d_ULA_nor,Nr_UCA_index,Nr_ULA_index,Nr_UCA,Nr_ULA,L,M,Nt);
        dev_h_delay=[dev_h_delay; transpose(Br_delay)];

        Br_angle_Phi = SEMI_spec_chan_derive_angle_Phi_UCyA(fading,delay,DOA_Phi,DOA_Theta,R_nor,d_ULA_nor,Nr_UCA_index,Nr_ULA_index,Nr_UCA,Nr_ULA,L,M,Nt);
        dev_h_angle_Phi=[dev_h_angle_Phi; transpose(Br_angle_Phi)];

        Br_angle_Theta = SEMI_spec_chan_derive_angle_Theta_UCyA(fading,delay,DOA_Phi,DOA_Theta,R_nor,d_ULA_nor,Nr_UCA_index,Nr_ULA_index,Nr_UCA,Nr_ULA,L,M,Nt);
        dev_h_angle_Theta=[dev_h_angle_Theta; transpose(Br_angle_Theta)];
    end
end

%% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters
G = [dev_h_fading,dev_h_delay,dev_h_angle_Phi,dev_h_angle_Theta]; 
%% ------------------------------------------------------------------------
 
X = [];
for ii = 1 : Nt
    X  = [X diag(ZC(:,ii))*FL];
end

[H, h_true] = gen_chan_specular(fading,delay,DOA_Phi,DOA_Theta,R_nor,d_ULA_nor,Nr_UCA,Nr_ULA,L,Nt);


%LAMBDA
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


% Partial derivative of LAMBDA w.r.t. h_i
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

N_total=4;
N_pilot=1;
N_data=N_total-N_pilot;
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
    Iop      = X_nga'*X_nga / sigmav2;
    Iop_full = N_total * Iop;
    CRB_op(snr_i) = abs(trace(pinv(Iop_full)));
%============================================
    %Only Pilot Specular   
    Iop_spec = G*G'*Iop_full*G*G';
    CRB_op_spec(snr_i) = abs(trace(pinv(Iop_spec)));
%============================================
%SemiBlind
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
%1 pilot
    I_SB_1               = N_data*I_D+N_pilot*Iop;
    CRB_SB_i_1           = pinv(I_SB_1);
    CRB_SB_1(snr_i)      = abs(trace(CRB_SB_i_1));
    
%2 pilots
    I_SB_2               = (N_data-1)*I_D+(N_pilot+1)*Iop;
    CRB_SB_i_2           = pinv(I_SB_2);
    CRB_SB_2(snr_i)      = abs(trace(CRB_SB_i_2));
    
%3 pilots
    I_SB_3               = (N_data-2)*I_D+(N_pilot+2)*Iop;
    CRB_SB_i_3           = pinv(I_SB_3);
    CRB_SB_3(snr_i)      = abs(trace(CRB_SB_i_3));
    
    clear I_D;
    toc
    
%============================================
%Semiblind Specular
%1 pilot
   I_SB_spec_1=G*G'*I_SB_1*G*G';
   CRB_SB_spec_i_1           = pinv(I_SB_spec_1);
   CRB_SB_spec_1(snr_i)      = abs(trace(CRB_SB_spec_i_1));  
   
%2 pilots
   I_SB_spec_2=G*G'*I_SB_2*G*G';
   CRB_SB_spec_i_2           = pinv(I_SB_spec_2);
   CRB_SB_spec_2(snr_i)      = abs(trace(CRB_SB_spec_i_2));  

%3 pilots
   I_SB_spec_3=G*G'*I_SB_3*G*G';
   CRB_SB_spec_i_3           = pinv(I_SB_spec_3);
   CRB_SB_spec_3(snr_i)      = abs(trace(CRB_SB_spec_i_3));    
   
end

%figure
semilogy(SNR,CRB_op,'-b')
hold on; semilogy(SNR,CRB_op_spec,'-r')

semilogy(SNR,CRB_SB_1,'-g')
hold on; semilogy(SNR,CRB_SB_spec_1,'-m')

semilogy(SNR,CRB_SB_2,'-b>')
hold on; semilogy(SNR,CRB_SB_spec_2,'-r+')


semilogy(SNR,CRB_SB_3,'-g>')
hold on; semilogy(SNR,CRB_SB_spec_3,'-m+')


grid on
ylabel('Normalized CRB')
xlabel('SNR(dB)')
legend('normal OP','spec OP','normal SB 1p','spec SB 1p','normal SB 2p','spec SB 2p','normal SB 3p','spec SB 3p')
%legend('normal_OP','spec_OP','normal_SB','spec_SB')
title(' ')
%axis([-10 20 1e-4 1e2])
hold on;