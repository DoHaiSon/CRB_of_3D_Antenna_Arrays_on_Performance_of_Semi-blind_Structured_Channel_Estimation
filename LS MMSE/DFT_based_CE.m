% channel_estimation.m
% for LS/DFT Channel Estimation with linear/spline interpolation
clear all; close all; clf

Nfft=64;
Ng=Nfft/8;
Nofdm=Nfft+Ng;
Nsym=128;
Nps=4;
Np=Nfft/Nps; % Pilot spacing and number of pilots per OFDM symbol
Nbps=4;
M=2^Nbps; % Number of bits per (modulated) symbol
mod_object = modem.qammod('M', M, 'SymbolOrder', 'gray');
demod_object = modem.qamdemod('M', M, 'SymbolOrder','gray');
Es=1;
A=sqrt(3/2/(M-1)*Es); % Signal energy and QAM normalization factor
SNR = 30; % dB
sq2=sqrt(2);
MSE = zeros(1, 6);
nose = 0;
for nsym=1:Nsym
    Xp = 2*(randn(1,Np)>0)-1; % Pilot sequence generation
    msgint=randint(1,Nfft-Np,M); % bit generation
    Data = A*modulate(mod_object, msgint);
    ip = 0;
    pilot_loc = [];
    for k=1:Nfft
        if mod(k,Nps)==1
            X(k)=Xp(floor(k/Nps)+1);
            pilot_loc=[pilot_loc k];
            ip = ip+1;
        else
            X(k) = Data(k-ip);
        end
    end
    x = ifft(X,Nfft);
    xt = [x(Nfft-Ng+1:Nfft) x]; % IFFT and add CP
    h = [(randn+1i*randn) (randn+1i*randn)/2]; % A (2-tap) channel
    H = fft(h,Nfft);
    ch_length=length(h); % True channel and its length
    H_power_dB = 10*log10(abs(H.*conj(H))); % True channel power in dB
    y_channel = conv(xt,h); % Channel path (convolution)
    yt = awgn(y_channel, SNR , 'measured');
    y = yt(Ng+1:Nofdm); % Remove CP
    Y = fft(y); % FFT
    for m=1:3
        if m==1
            H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'linear');
            method='LS-linear'; % LS estimation with linear interpolation
        elseif m==2
            H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'spline');
            method='LS-spline'; % LS estimation with spline interpolation
        else
            H_est = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR);
            method='MMSE'; % MMSE estimation
        end
        H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
        h_est = ifft(H_est);
        h_DFT = h_est(1:ch_length);
        H_DFT = fft(h_DFT,Nfft); % DFT-based channel estimation
        H_DFT_power_dB = 10*log10(abs(H_DFT.*conj(H_DFT)));
        if nsym==1
            subplot(319+2*m);
            plot(H_power_dB,'b');
            hold on;
            plot(H_est_power_dB,'r:+');
            legend('True Channel',method);
            
            subplot(320+2*m);
            plot(H_power_dB,'b');
            hold on;
            plot(H_DFT_power_dB,'r:+');
            legend('True Channel',[method ' with DFT']);
        end
        MSE(m) = MSE(m) + (H-H_est)*(H-H_est)';
        MSE(m+3) = MSE(m+3) + (H-H_DFT)*(H-H_DFT)';
    end
    Y_eq = Y./H_est;
    ip = 0;
    for k=1:Nfft
        if mod(k,Nps)==1
            ip=ip+1;
        else
            Data_extracted(k-ip)=Y_eq(k);
        end
    end
    msg_detected = demodulate(demod_object,Data_extracted/A);
    nose = nose + sum(msg_detected~=msgint);
    MSEs = MSE/(Nfft*Nsym);
end