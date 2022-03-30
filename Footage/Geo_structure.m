clc
close all;
clear all;

Nr_UCA  = 8;         % number of receive antennas of UCA
Nr_ULA  = 1;         % number of receive antennas of ULA

d_ULA_nor   = 0.5;
d_UCA_nor   = 0.5;
R_nor       = 0.5 * d_UCA_nor/sin(pi/Nr_UCA);

ULA_nor                 = zeros(Nr_ULA, Nr_UCA);
position_elements_nor   = zeros(3, Nr_ULA, Nr_UCA);

for Nr_ULA_index=1:Nr_ULA
    for Nr_UCA_index=1:Nr_UCA
        if Nr_ULA == 1
            position_elements_nor(1, Nr_ULA_index, Nr_UCA_index) = (Nr_UCA_index-1) * d_UCA_nor;
            position_elements_nor(2, Nr_ULA_index, Nr_UCA_index) = 0;
            position_elements_nor(3, Nr_ULA_index, Nr_UCA_index) = 0;
        else
            position_elements_nor(1, Nr_ULA_index, Nr_UCA_index) = R_nor * sin((Nr_UCA_index-1)*(2*pi/Nr_UCA)) ;         % x
            position_elements_nor(2, Nr_ULA_index, Nr_UCA_index) = R_nor * cos((Nr_UCA_index-1)*(2*pi/Nr_UCA)) ;         % y
            position_elements_nor(3, Nr_ULA_index, Nr_UCA_index) = (Nr_ULA_index-1) * d_ULA_nor;                         % z
        end
%         ULA_nor(Nr_ULA_index, Nr_UCA_index) = (Nr_UCA_index-1)*2*pi/Nr_UCA + rot_nor * (Nr_ULA_index-1)*2*pi/Nr_UCA;    % polar axis
    end
end

x = squeeze(position_elements_nor(1, :, :));
% x(:, end + 1) = x(:, 1);
y = squeeze(position_elements_nor(2, :, :));
% y(:, end + 1) = y(:, 1);
z = squeeze(position_elements_nor(3, :, :));
% z(:, end + 1) = z(:, 1);

x = x.';
y = y.';
z = z.';

plot3(x(1,:), y(1, :), z(1, :));
hold on;
% plot3(x(2,:), y(2, :), z(2, :));
% hold on;
% plot3(x(3,:), y(3, :), z(3, :));
% hold on;
% plot3(x(4,:), y(4, :), z(4, :));
grid on