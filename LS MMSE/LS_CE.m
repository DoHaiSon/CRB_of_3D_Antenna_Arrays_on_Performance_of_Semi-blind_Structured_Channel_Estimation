function [LS_est] = LS_CE(Y, Xp, Np)
% LS channel estimation function
% Inputs:
% Y = Frequency-domain received signal
% Xp = Pilot signal
% Np = Number of pilots
% output:
% LS_est = LS Channel estimate
    k=1:Np;
    LS_est(k) = Y(k)./Xp(k); % LS channel estimation
end