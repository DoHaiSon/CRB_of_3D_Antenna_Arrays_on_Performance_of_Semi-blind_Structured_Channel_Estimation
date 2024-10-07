clear all;
close all;
clc; 

addpath('./funcs');

%% Declear

Nt      = 8;                    % number of transmit antennas

Nr_UCA  = 16;                   % UCA
Nr_ULA  = 4;                    % UCyA
Nr      = Nr_UCA * Nr_ULA;      % ULA

d_H     = 1;                    % Number of propagation paths

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

H_f = [];
H_tmp_f = [];

for t = 1: 30
    H = sqrt(1/2) * (normrnd(0, 1, Nr, Nt) + 1i*normrnd(0, 1, Nr, Nt));
    H_f = [H_f, H(:)];
    %% Channel generation
    gamma       = zeros(d_H, Nt);      % complex gain
    ZOA         = zeros(d_H, Nt);      % ffset ZOA ray
    AOA         = zeros(d_H, Nt);      % offset AOA ray
    
    for nt = 1 : Nt
%         gamma(:, nt)   = sqrt(0.5)*(normrnd(0,1,1,d_H) + 1j*normrnd(0,1,1,d_H));
%         while min(abs(gamma(:, nt)))<0.6
%               gamma(:, nt)   = sqrt(0.5)*(normrnd(0,1,1,d_H) + 1j*normrnd(0,1,1,d_H));
%         end
        
        AOA(:, nt)     = random('unif',0,1,1,d_H)*pi-pi/2; 
%         dist=pdist(vec(AOA(:, nt)),'euclid');
%         while min(dist) < pi/18
%               AOA(:, nt)=random('unif',0,1,1,d_H)*pi-pi/2;
%               dist=pdist(vec(AOA(:, nt)),'euclid');
%         end
    
        ZOA(:, nt)     = random('unif',0,1,1,d_H)*pi-pi/2; 
%         dist=pdist(vec(ZOA(:, nt)),'euclid');
%         while min(dist) < pi/18
%               ZOA(:, nt)=random('unif',0,1,1,d_H)*pi-pi/2;
%               dist=pdist(vec(ZOA(:, nt)),'euclid');
%         end
    end 


    r = 0;
    H_tmp = zeros(Nr, Nt);
    for Nr_ULA_index=1:Nr_ULA   
        for Nr_UCA_index=1:Nr_UCA
            r=r+1;
            tmp = zeros(d_H, Nt);
            for jj = 1 : Nt
                    r_x = sin(ZOA(1,jj)) * cos(AOA(1,jj));
                    r_y = sin(ZOA(1,jj)) * sin(AOA(1,jj));
                    r_z = cos(ZOA(1,jj));
                    H_tmp(r, jj) = H(r, jj) / exp(-1i*2*pi* ...
                         (UCyA_elements_nor(1, Nr_ULA_index, Nr_UCA_index)*r_x ...
                        + UCyA_elements_nor(2, Nr_ULA_index, Nr_UCA_index)*r_y ...
                        + UCyA_elements_nor(3, Nr_ULA_index, Nr_UCA_index)*r_z));
            end
        end
    end
    H_tmp_f = [H_tmp_f, H_tmp(:)];
end
%     [H, h_vec] = gen_chan_specular_UCyA( gamma, AOA, ZOA, UCyA_elements_nor, Nt, Nr_UCA, Nr_ULA, d_H );

    %% Generate Beta matrix (MT x T)
%     gamma_vec = gamma(:);
%     Gamma = repmat(gamma_vec, 1, Nt );
%     tmp = zeros(d_H * Nt, Nt);
%     for nt=1:Nt
%         for j= (nt-1)*d_H + 1 : nt*d_H
%             tmp(j, nt) = Gamma(j, nt);
%         end
%     end
%     Gamma = tmp;

    %% Generate Steering matrix (L x MT)
%     r = 0;
%     Steering = zeros(Nr, Nt * d_H);
%     for Nr_ULA_index=1:Nr_ULA   
%         for Nr_UCA_index=1:Nr_UCA
%             r=r+1;
%             tmp = zeros(d_H, Nt);
%             for jj = 1 : Nt
%                 for nn = 1 : d_H 
%                     r_x = sin(ZOA(nn,jj)) * cos(AOA(nn,jj));
%                     r_y = sin(ZOA(nn,jj)) * sin(AOA(nn,jj));
%                     r_z = cos(ZOA(nn,jj));
%                     tmp(nn, jj) = exp(-1i*2*pi* ...
%                          (UCyA_elements_nor(1, Nr_ULA_index, Nr_UCA_index)*r_x ...
%                         + UCyA_elements_nor(2, Nr_ULA_index, Nr_UCA_index)*r_y ...
%                         + UCyA_elements_nor(3, Nr_ULA_index, Nr_UCA_index)*r_z));
%                 end
%             end
%             Steering(r, :) = tmp(:);
%         end
%     end
% 
%     H_1 = Steering * Gamma;

%     [H, h_vec] = gen_chan_specular_ULA(gamma, AOA, ZOA, ULA_elements_nor,  Nt, Nr, d_H);
%     H_vec = [H_vec; h_vec];
% end

H_vec = H_f(:);

H_vec_r = real(H_vec);
H_vec_i = imag(H_vec);

histfit(H_vec_r);

pd_r = fitdist(H_vec_r,'Normal');
mu = pd_r.mu;
sigma = pd_r.sigma;
legend({'Histogram of $h_{l, t}$', strcat('$\mathcal{N}(', num2str(mu), ',', num2str(sigma), ')$')}, 'Interpreter','latex')
xlim([-3 3]);
grid minor;
set(gcf,'color','w');
ax = get(gca,'XTickLabel');
set(gca,'FontName','Times','fontsize',12);
ylabel('Frequency', 'FontSize', 16, 'Interpreter','latex');
xlabel('$\Re(h_{l, t})$', 'FontSize', 16, 'Interpreter','latex');

figure();
H_vec = H_tmp_f(:) / 1.1;

H_vec_r = real(H_vec);
H_vec_i = imag(H_vec);

histfit(H_vec_r);

pd_r = fitdist(H_vec_r,'Normal');
mu_1 = pd_r.mu;
sigma_1 = pd_r.sigma;

legend({'Histogram of $\hat{\beta_{l, t}}$', strcat('$\mathcal{N}(', num2str(mu_1), ',', num2str(sigma_1), ')$')}, 'Interpreter','latex')
xlim([-3 3]);
grid minor;
set(gcf,'color','w');
ax = get(gca,'XTickLabel');
set(gca,'FontName','Times','fontsize',12);
ylabel('Frequency', 'FontSize', 16, 'Interpreter','latex');
xlabel('$\Re(\hat{\beta_{l, t}})$', 'FontSize', 16, 'Interpreter','latex');