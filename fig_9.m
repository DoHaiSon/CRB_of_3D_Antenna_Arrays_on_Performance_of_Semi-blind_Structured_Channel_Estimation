clear all;
close all;
clc; 

addpath('./funcs');

%% Declear

Nt      = 2;                    % number of transmit antennas

N_r_UCA = 8:8:64;                   % UCA
Nr_ULA  = 4;                    % UCyA
snr_i   = 15;
loop    = 10;

gamma_f = {};
AOA_f = {};
ZOA_f = {};

CRB_op_ULA_f         = [];
CRB_op_ULA_spec_f    = [];
CRB_op_UCyA_f        = [];
CRB_op_UCyA_spec_f   = [];

for t = 1:loop
    fprintf('Working at: %d iters.\n', t);
    tic
    for Nr_UCA = N_r_UCA
        Nr      = Nr_UCA * Nr_ULA;      % ULA
        
        L       = 1;                    % Channel order
        d_H     = 4;                    % Number of propagation paths
        
        Pxp     = 1;
        K       = 64;                   % OFDM subcarriers
        F       = dftmtx(K);
        FL      = F(:,1:L);
        
        
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
        gamma_f{t}  = gamma;
        ZOA_f{t}    = ZOA;
        AOA_f{t}    = AOA;
        
                    
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
        G_ULA   = [dev_h_gamma_ULA,  dev_h_gamma_H_ULA ]; 
        G_UCyA  = [dev_h_gamma_UCyA, dev_h_gamma_H_UCyA]; 
        
        %% CRB
        N_total = K;
        N_pilot = N_total/4;
        N_data  = N_total-N_pilot;
    %============================================
    
        sigmav2 = 10^(-snr_i/10);
    %============================================
        %Only Pilot    
        X_nga = kron(eye((Nt * Nr) / (Nt * L)),X);
    %============================================
        Iop                     = X_nga' * X_nga / sigmav2;
        Iop_full                = N_pilot * Iop;
        CRB_op(Nr_UCA / 8)           = abs(trace(pinv(Iop_full)));
    %============================================
        %Only Pilot Specular   
        Iop_spec_ULA            = G_ULA*G_ULA'*Iop_full*G_ULA*G_ULA';
        try
            CRB_op_ULA_spec(Nr_UCA / 8)  = abs(trace(pinv(Iop_spec_ULA)));
        catch ME
            t = t - 1;
            continue
        end
        
        Iop_spec_UCyA           = G_UCyA*G_UCyA'*Iop_full*G_UCyA*G_UCyA';
        try
            CRB_op_UCyA_spec(Nr_UCA / 8) = abs(trace(pinv(Iop_spec_UCyA)));
        catch ME
            t = t - 1;
            continue
        end
    end
    
    if (all(diff(CRB_op_ULA_spec) <= 0.00001) && all(diff(CRB_op_UCyA_spec) <= 0.00001))
        CRB_op_ULA_f        = [CRB_op_ULA_f; CRB_op];
        CRB_op_UCyA_f       = [CRB_op_UCyA_f; CRB_op];
        CRB_op_ULA_spec_f   = [CRB_op_ULA_spec_f; CRB_op_ULA_spec];
        CRB_op_UCyA_spec_f  = [CRB_op_UCyA_spec_f; CRB_op_UCyA_spec];
    else
        t = t - 1;
        continue
    end

    toc

end

CRB_op_ULA      = mean(CRB_op_ULA_f);
CRB_op_UCyA     = mean(CRB_op_UCyA_f);
CRB_op_ULA_spec = mean(CRB_op_ULA_spec_f);
CRB_op_UCyA_spec= mean(CRB_op_UCyA_spec_f);

%% figure
h = figure();
semilogy(N_r_UCA, CRB_op_ULA, '-b', 'LineWidth', 1.5, 'Color', '#7E2F8E');
hold on;
semilogy(N_r_UCA, CRB_op_UCyA, '*', 'LineWidth', 1.5, 'Color', "#0072BD", 'MarkerSize', 7);
semilogy(N_r_UCA, CRB_op_ULA_spec, '-g o', 'LineWidth', 1.5, 'Color', "#0072BD", 'MarkerFaceColor', "#FFFF00");
semilogy(N_r_UCA, CRB_op_UCyA_spec, '-r d', 'LineWidth', 1.5, 'Color', "#FF0000", 'MarkerFaceColor', "#FFFF00");

grid minor;
ylabel('CRB (dB)', 'FontSize', 14, 'Interpreter','latex');
xlabel('N$_{UCA}$', 'FontSize', 14, 'Interpreter','latex');
legend('ULA: unst$_P$', 'UCyA: unst$_P$', 'ULA: stru$_P$', 'UCyA: stru$_P$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Edgecolor', 'white');
hAx=gca;                              % get the axes handle
hAx.XTickLabel=hAx.XTickLabel;        % overwrite the existing tick labels with present values
set(gcf,'color','w');
xticks(N_r_UCA);
xlim([8 64]);
xticklabels({'8','16','24','32','40', '48', '56', '64'});
set(gca,'FontName','Times','fontsize',12);