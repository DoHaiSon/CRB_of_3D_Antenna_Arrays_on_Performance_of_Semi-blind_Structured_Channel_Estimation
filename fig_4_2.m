clear all;
close all;
clc; 

addpath('./funcs');

%% Declear

Nt      = 2;                    % number of transmit antennas

Nr_UCA  = 24;                   % UCA
N_r_ULA = 1:10;                 % UCyA
snr_i   = 5;
loop    = 100;
loop_i  = 1;

CRB_SB_ULA_f         = [];
CRB_SB_ULA_spec_f    = [];
CRB_SB_UCyA_spec_f   = [];

I_D_f = {};

for Nr_ULA = N_r_ULA
    I_D_f{Nr_ULA} = gather(getfield(load('./I_D/SNR_5.mat'), 'I_D'));
end

while true

    fprintf('Working at iter: %d.\n', loop_i);
    tic
    for Nr_ULA = N_r_ULA

        Nr      = Nr_UCA * Nr_ULA;      % ULA
        
        L       = 1;                    % Channel order
        d_H     = 4;                    % Number of propagation paths
        
        Pxp     = 1;
        K       = 64;                   % OFDM subcarriers
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
        U = 1:2:100;
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
        gamma       = zeros(d_H, Nt);      % complex gain
        ZOA         = zeros(d_H, Nt);      % ffset ZOA ray
        AOA         = zeros(d_H, Nt);      % offset AOA ray
        
        for nt = 1 : Nt
            gamma(:, nt)   = sqrt(0.5)*(normrnd(0,1,1,d_H) + 1j*normrnd(0,1,1,d_H));
            while min(abs(gamma(:, nt)))<0.6
                  gamma(:, nt)   = sqrt(0.5)*(normrnd(0,1,1,d_H) + 1j*normrnd(0,1,1,d_H));
            end
            
            AOA(:, nt)     = random('unif',0,1,1,d_H)*pi-pi/2; 
            dist=pdist(vec(AOA(:, nt)),'euclid');
            while min(dist) < pi/18
                  AOA(:, nt)=random('unif',0,1,1,d_H)*pi-pi/2;
                  dist=pdist(vec(AOA(:, nt)),'euclid');
            end
    
            ZOA(:, nt)     = random('unif',0,1,1,d_H)*pi-pi/2; 
            dist=pdist(vec(ZOA(:, nt)),'euclid');
            while min(dist) < pi/18
                  ZOA(:, nt)=random('unif',0,1,1,d_H)*pi-pi/2;
                  dist=pdist(vec(ZOA(:, nt)),'euclid');
            end
        end     
        
                    
        %% Derivative
        dev_h_gamma_ULA     = [];
        dev_h_gamma_H_ULA   = [];
        dev_h_ZOA_ULA       = [];
        dev_h_AOA_ULA       = [];
        
        dev_h_gamma_UCyA    = [];
        dev_h_gamma_H_UCyA  = [];
        dev_h_ZOA_UCyA      = [];
        dev_h_AOA_UCyA      = [];
        
        for Nr_index=1:Nr
            Br_gamma                = SEMI_spec_chan_derive_gamma(   gamma, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), d_H, Nt);
            dev_h_gamma_ULA         = [dev_h_gamma_ULA; transpose(Br_gamma)];
    
            Br_gamma_H              = SEMI_spec_chan_derive_gamma_H( gamma, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), d_H, Nt);
            dev_h_gamma_H_ULA       = [dev_h_gamma_H_ULA; transpose(Br_gamma_H)];
            
            Br_angle_ZOA            = SEMI_spec_chan_derive_ZOA(     gamma, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), d_H, Nt);
            dev_h_ZOA_ULA           = [dev_h_ZOA_ULA; transpose(Br_angle_ZOA)];
    
            Br_angle_AOA            = SEMI_spec_chan_derive_AOA(     gamma, AOA, ZOA, ULA_elements_nor(:, 1, Nr_index), d_H, Nt);
            dev_h_AOA_ULA           = [dev_h_AOA_ULA; transpose(Br_angle_AOA)];
        end
        
        for Nr_ULA_index=1:Nr_ULA
            for Nr_UCA_index=1:Nr_UCA
                Br_gamma            = SEMI_spec_chan_derive_gamma(   gamma, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), d_H, Nt);
                dev_h_gamma_UCyA    = [dev_h_gamma_UCyA; transpose(Br_gamma)];
    
                Br_gamma_H          = SEMI_spec_chan_derive_gamma_H( gamma, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), d_H, Nt);
                dev_h_gamma_H_UCyA  = [dev_h_gamma_H_UCyA; transpose(Br_gamma_H)];
    
                Br_angle_ZOA        = SEMI_spec_chan_derive_ZOA(     gamma, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), d_H, Nt);
                dev_h_ZOA_UCyA      = [dev_h_ZOA_UCyA; transpose(Br_angle_ZOA)];
    
                Br_angle_AOA        = SEMI_spec_chan_derive_AOA(     gamma, AOA, ZOA, UCyA_elements_nor(:, Nr_ULA_index, Nr_UCA_index), d_H, Nt);
                dev_h_AOA_UCyA      = [dev_h_AOA_UCyA; transpose(Br_angle_AOA)];
            end
        end
    
        %% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters
        G_ULA   = [dev_h_gamma_ULA,  dev_h_gamma_H_ULA,  dev_h_ZOA_ULA,  dev_h_AOA_ULA]; 
        G_UCyA  = [dev_h_gamma_UCyA, dev_h_gamma_H_UCyA, dev_h_ZOA_UCyA, dev_h_AOA_UCyA]; 
    
        %% Partial lambda
        [H_ULA,  h_true_ULA,  LAMBDA_ULA,  partial_LAMBDA_ULA]  = partial_LAMBDA_ULA_f ( gamma, AOA, ZOA, ULA_elements_nor,  Nt, Nr, d_H, FL);
        
        [H_UCyA, h_true_UCyA, LAMBDA_UCyA, partial_LAMBDA_UCyA] = partial_LAMBDA_UCyA_f( gamma, AOA, ZOA, UCyA_elements_nor, Nt, Nr_UCA, Nr_ULA, d_H, FL );
        
        %% CRB
        N_total = K;
        N_pilot = K/4;
        N_data  = N_total-N_pilot;
    %============================================

    %% Only Pilot    
        sigmav2 = 10^(-snr_i/10);
        X_nga = kron(eye(Nr),X);
        Iop   = X_nga' * X_nga / sigmav2;
        
    %% SemiBlind
        I_D      = I_D_f{Nr_ULA};

        I_D = triu(repmat(I_D, Nr*Nt, Nr*Nt));
    %============================================
        %Semiblind Normal
        I_SB               = N_data*I_D+N_pilot*Iop;
        CRB_SB_i           = pinv(I_SB);
        CRB_SB_ULA(Nr_ULA)  = abs(trace(CRB_SB_i));

    %============================================
        %Semiblind Specular
        I_SB_spec               = G_ULA*G_ULA'*I_SB*G_ULA*G_ULA';
        CRB_SB_spec_i           = pinv(I_SB_spec);
        CRB_SB_ULA_spec(Nr_ULA)  = abs(trace(CRB_SB_spec_i));   
  
        
 %% ========================================================================================
        clear I_D;    

        %Semiblind Specular
        I_SB_spec               = G_UCyA*G_UCyA'*I_SB*G_UCyA*G_UCyA';
        CRB_SB_spec_i           = pinv(I_SB_spec);
        CRB_SB_UCyA_spec(Nr_ULA) = abs(trace(CRB_SB_spec_i));
    end

    if (sum(isnan(CRB_SB_ULA_spec)) == 0 && sum(isnan(CRB_SB_UCyA_spec)) == 0 ...
            && all(diff(CRB_SB_ULA_spec) <= 0.000001) ...
            && all(diff(CRB_SB_UCyA_spec) <= 0.000001))
        CRB_SB_ULA_f        = [CRB_SB_ULA_f; CRB_SB_ULA];
        CRB_SB_ULA_spec_f   = [CRB_SB_ULA_spec_f; CRB_SB_ULA_spec];
        CRB_SB_UCyA_spec_f  = [CRB_SB_UCyA_spec_f; CRB_SB_UCyA_spec];
    else
        continue
    end

    loop_i = loop_i + 1;
    toc
    if(loop_i > loop)
        break
    end
end

CRB_SB_ULA      = mean(CRB_SB_ULA_f);
CRB_SB_ULA_spec = mean(CRB_SB_ULA_spec_f);
CRB_SB_UCyA_spec= mean(CRB_SB_UCyA_spec_f);

%% figure
h = figure();
semilogy(N_r_ULA, CRB_SB_ULA, '-b', 'LineWidth', 1.5, 'Color', '#7E2F8E');
hold on;
semilogy(N_r_ULA, CRB_SB_ULA_spec, '-g o', 'LineWidth', 1.5, 'Color', "#0072BD", 'MarkerFaceColor', "#FFFF00");
semilogy(N_r_ULA, CRB_SB_UCyA_spec, '-r d', 'LineWidth', 1.5, 'Color', "#FF0000", 'MarkerFaceColor', "#FFFF00");

grid minor;
ylabel('CRB (dB)', 'FontSize', 14, 'Interpreter','latex');
xlabel('N$_{3D}$', 'FontSize', 14, 'Interpreter','latex');
legend('unst$_{SB}$', 'ULA: stru$_{SB}$', 'UCyA: stru$_{SB}$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Edgecolor', 'white');
hAx=gca;                              % get the axes handle
hAx.XTickLabel=hAx.XTickLabel;        % overwrite the existing tick labels with present values
set(gcf,'color','w');
xticks(N_r_ULA);
xticklabels({'1','2','3','4','5', '6', '7', '8', '9', '10'});
set(gca,'FontName','Times','fontsize',12);