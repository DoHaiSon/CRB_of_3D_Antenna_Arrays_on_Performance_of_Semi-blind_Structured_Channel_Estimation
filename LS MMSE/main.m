% channel_estimation.m
% for LS/DFT Channel Estimation with linear/spline interpolation
clear all; close all;

Nfft            = 48;               % Occ carriers
Ng              = Nfft/4;           % CP len
Nofdm           = Nfft+Ng;          % OFDM symbol len
Nbit            = Nfft;             % Seq
Nps             = 1;                % Pilot spacing
Np              = Nfft/Nps;         % Number of pilots per OFDM symbol
Nbps            = 1;                % Number of bits per symbol
M               = 2^Nbps;           % Number of waveforms

MSE_LS_f = [];
MSE_LS_DFT_f = [];
ser_f = [];
for t = 1:1000
    %% Pilot sequence generation
    Xp      = (2*(randn(1, Np)>0)-1) * (0.7 + 0.7i); 
    %% Data Bit generation
    msgint  = randint(1, Nbit, M);        
    x_d     = qammod(msgint, M);

    x_p     = ifft(Xp,Nfft);
    xt_p    = [x_p(Nfft-Ng+1:Nfft) x_p]; % IFFT and add CP
    x_d     = ifft(x_d,Nfft);
    xt_d    = [x_d(Nfft-Ng+1:Nfft) x_d]; % IFFT and add CP

    h       = 1/sqrt(2) * (rand(1, 3) + 1i*rand(1, 3)); % A (3-tap) channel
    H       = fft(h,Nfft);
    ch_length=length(h); % True channel and its length
    H_power_dB = 10*log10(abs(H.*conj(H))); % True channel power in dB

    %% Channel path (convolution)
    y_channel       = conv(xt_p, h); 
    y_d_channel     = conv(xt_d, h);

    %% SNR
    SNR = -20:1:20;
    MSE_LS = [];
    MSE_LS_DFT = [];
    ser = [];

    for snr = SNR
        yt_p      = awgn(y_channel, snr, 'measured');
        yt_d      = awgn(y_d_channel, snr, 'measured');

        %% Remove CP
        y_p         = yt_p(Ng+1:Nofdm);
        Y_p         = fft(y_p); % FFT

        y_d         = yt_d(Ng+1:Nbit + Ng);
        Y_d         = fft(y_d); % FFT

        H_est     = LS_CE(Y_p, Xp, Np);
        method    = 'LS'; % LS estimation 
        H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
        h_est     = ifft(H_est);
        h_DFT     = h_est(1:ch_length);
        H_DFT     = fft(h_DFT,Nfft); % DFT-based channel estimation
        H_DFT_power_dB = 10*log10(abs(H_DFT.*conj(H_DFT)));

        %% Pilot H_est
    %     figure();
    %     subplot(1,2,1);
    %     plot(H_power_dB,'b');
    %     hold on;
    %     plot(H_est_power_dB,'r:+');
    %     legend('True Channel',method);
    % 
    %     subplot(1,2,2);
    %     plot(H_power_dB,'b');
    %     hold on;
    %     plot(H_DFT_power_dB,'r:+');
    %     legend('True Channel',[method ' with DFT']);

        %% MSE H
        MSE_LS      = [MSE_LS, (H-H_est)*(H-H_est)'];
        MSE_LS_DFT  = [MSE_LS_DFT, (H-H_DFT)*(H-H_DFT)'];

        x_est   = [];
        for i = 1:Nbit / Nfft
            x_est = [x_est, Y_d((i-1)*Nfft + 1:i*Nfft) ./ H_est(1:Nfft)];
        end

        msg_detected = qamdemod(x_est, M);

        err = sum(msgint ~= msg_detected) / Nbit * 100;
        ser = [ser, err];
    end
    MSE_LS_f = [MSE_LS_f; MSE_LS];
    MSE_LS_DFT_f = [MSE_LS_DFT_f; MSE_LS_DFT];
    ser_f = [ser_f; ser];
end

MSE_LS = mean(MSE_LS_f);
MSE_LS_DFT = mean(MSE_LS_DFT_f);
ser  = mean(ser_f);

figure();
semilogy(SNR, MSE_LS);
hold on;
semilogy(SNR, MSE_LS_DFT);
grid on;
ylabel('MSE h');
xlabel('SNR(dB)');
legend('MSE LS', 'MSE LS DFT');

figure();
semilogy(ser);
grid on;
ylabel('Symbol error rate (SER)');
xlabel('SNR(dB)');