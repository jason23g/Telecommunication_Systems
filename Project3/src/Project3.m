clear all;
close all;

M = 16;%16-QAM
N = 200;
A = 1;
A_pulse = 4;
T = 0.01;
over = 10;
Ts = T/over ;
DT = Ts;
Fs = 1/Ts;
a = 0.5;
Nf = 2048;
F0 = 200;
K = 1000;
SNR_dB = [0:2:16];
F = [-Fs/2:Fs/Nf:Fs/2 - Fs/Nf];

Tot_E_symbol = zeros(1,length(SNR_dB));
Tot_E_bits = zeros(1,length(SNR_dB));

Theoretical_Tot_E_symbol = zeros(1,length(SNR_dB));
Theoretical_Tot_E_bits = zeros(1,length(SNR_dB));

for j = 1:1:length(SNR_dB)

Average_E_symbols  = 0;
Average_E_bits = 0;

for i = 1:1:K
%----------------------------PART A--------------------%
%%%%%%%%%%%%%%%%%%Part 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = (sign(randn(4*N,1)) + 1)/2;%Generating a sequence of bits

%%%%%%%%%%%%%%%%%%Part 2 and 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Seperating the sequence b of bits in two sequences of bits
b1 = b(1 :1 :2*N);

b2 =  b((2*N + 1) :1 :4*N);

%Creating sequences of symbols from 4-PAM alphabet
X_I = bits_to_4_PAM(b1,A);

X_Q = bits_to_4_PAM(b2,A);

%
X_I_d = (1/Ts)*upsample(X_I,over);

time_of_X_I_d = [0:Ts:(N*T - Ts)];

X_Q_d = (1/Ts)*upsample(X_Q,over);

time_of_X_Q_d = [0:Ts:(N*T - Ts)];

%%%%%%%%%%%%%%%%%%%Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creating a SRRC pulse
[phi, t] = srrc_pulse(T, Ts, A_pulse, a);

%Modulating the two sequences X_I_d and X_Q_d
X_i = Ts*conv(phi,X_I_d);

time_conv_X_i = [time_of_X_I_d(1) + t(1): DT : time_of_X_I_d(end) + t(end)];

X_q = Ts*conv(phi,X_Q_d);

time_conv_X_q = [time_of_X_Q_d(1) + t(1): DT : time_of_X_Q_d(end) + t(end)];

%Plotting X_I and X_Q
% figure(1);
% subplot(2,1,1);
% plot(time_conv_X_i,X_i);
% xlabel('Time(s)');
% ylabel('X_I');
% title('Convolution of X_I and srrc pulse \phi');
% 
% subplot(2,1,2);
% plot(time_conv_X_q,X_q);
% xlabel('Time(s)');
% ylabel('X_Q');
% title('Convolution of X_Q and srrc pulse \phi');

%Creating periodograms of X_i and X_q
t_total_X_i = time_conv_X_i(end)-time_conv_X_i(1);
P_X_i_F = (abs(fftshift(fft(X_i,Nf)*Ts)).^2)./(t_total_X_i);

t_total_X_q = time_conv_X_q(end)-time_conv_X_q(1);
P_X_q_F = (abs(fftshift(fft(X_q,Nf)*Ts)).^2)./(t_total_X_q);

%Plotting periodograms
% figure(2);
% subplot(2,1,1);
% plot(F,P_X_i_F);
% xlabel('Frequency(Hz)');
% ylabel('P{_{X_I}_F}');
% title('Periodogram of X_I');
% 
% subplot(2,1,2);
% plot(F,P_X_q_F);
% xlabel('Frequency(Hz)');
% ylabel('P{_{X_Q}_F}');
% title('Periodogram of X_Q');

%%%%%%%%%%%%%%%%%%%Part 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiplying X_i and X_q with the appropriate carrier 
X_I_mod = 2*X_i.*cos(2*pi*F0.*time_conv_X_i);

X_Q_mod = -2*X_q.*sin(2*pi*F0.*time_conv_X_q);

%Plotting X_I_mod and X_Q_mod
% figure(3);
% subplot(2,1,1);
% plot(time_conv_X_i,X_I_mod);
% xlabel('Time(s)');
% ylabel('{X_I}^{mod}');
% title('Modulation of X_I');
% 
% subplot(2,1,2);
% plot(time_conv_X_q,X_Q_mod);
% xlabel('Time(s)');
% ylabel('{X_Q}^{mod}');
% title('Modulation of X_Q');

%Creating periodograms of X_I_mod and X_Q_mod
P_X_i_mod_F = (abs(fftshift(fft(X_I_mod,Nf)*Ts)).^2)./(t_total_X_i);

P_X_q_mod_F = (abs(fftshift(fft(X_Q_mod,Nf)*Ts)).^2)./(t_total_X_q);
%Plotting periodograms of X_I_mod and X_Q_mod
% figure(4);
% subplot(2,1,1);
% plot(F,P_X_i_mod_F);
% xlabel('Frequency(Hz)');
% ylabel('P_{{X_I}^{mod}}_F');
% title('Periodogram of modulated X_I');
% 
% subplot(2,1,2);
% plot(F,P_X_q_mod_F);
% xlabel('Frequency(Hz)');
% ylabel('P_{{X_Q}^{mod}}_F');
% title('Periodogram of modulated X_Q');

%%%%%%%%%%%%%%%%%%%Part 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating X_mod which is the sum of X_I_mod and X_Q_mod
X_mod = X_I_mod + X_Q_mod ;

%Creating periodogram of X_mod
t_total_X = time_conv_X_i(end)-time_conv_X_i(1);
P_X_mod_F = (abs(fftshift(fft(X_mod,Nf)*Ts)).^2)./(t_total_X);

%Plotting X_mod and its periodogram
% figure(5);
% plot(time_conv_X_i,X_mod);
% xlabel('Time(s)');
% ylabel('X^{mod}');
% title('Modulated X');
% 
% figure(6);
% plot(F,P_X_mod_F);
% xlabel('Frequency(Hz)');
% ylabel('P_{X^{mod}}_F');
% title('Periodogram of modulated X');

%%%%%part 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = sqrt((10 * (A^2))/(Ts* (10^(SNR_dB(j)/10))));%Variance of Gaussian Noise

%Creating Noise
 Noise = sigma*randn(1,length(X_mod));

%Creating Y which is the sequence of symbols affected by Noise
Y = X_mod + Noise;

%%%%%part 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiplying Y with the appropriate carrier, creating Y_I_demod 
Y_I_demod = Y.*cos(2*pi*F0.*time_conv_X_i);

%Multiplying Y with the appropriate carrier, creating Y_Q_demod 
Y_Q_demod = Y.*(-sin(2*pi*F0.*time_conv_X_q));


%Plotting Y_i_demod and its periodogram
% figure(7);
% subplot(2,1,1);
% plot(time_conv_X_i,Y_I_demod);
% xlabel('Time(s)');
% ylabel('Y{_I}^{demod}');
% title('Demodulation of Y_I');
% 
% subplot(2,1,2);
% plot(time_conv_X_q,Y_Q_demod);
% xlabel('Time(s)');
% ylabel('Y{_Q}^{demod}');
% title('Demodulation of Y_Q');

%Creating periodogram of Y_i_demod
t_total_X_i = time_conv_X_i(end) -  time_conv_X_i(1);
P_Y_I_demod_F = (abs(fftshift(fft(Y_I_demod,Nf)*Ts)).^2)./(t_total_X_i);

%Creating periodogram of Y_q_demod
t_total_X_q = time_conv_X_q(end) -  time_conv_X_q(1);
P_Y_Q_demod_F = (abs(fftshift(fft(Y_Q_demod,Nf)*Ts)).^2)./(t_total_X_q);

%Plotting Y_q_demod and its periodogram
% figure(8);
% subplot(2,1,1);
% plot(F,P_Y_I_demod_F);
% xlabel('Frequency(Hz)');
% ylabel('P_{Y{_I}^{demod}}_F');
% title('Periodogram of demodulated X_I');
% 
% subplot(2,1,2);
% plot(F,P_Y_Q_demod_F);
% xlabel('Frequency(Hz)');
% ylabel('P_{Y{_Q}^{demod}}_F');
% title('Periodogram of demodulated X_Q');

%%%%%part 10%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Demodulating signal Y_I_mod
Y_I = Ts*conv(phi,Y_I_demod);
time_of_X_i = [time_conv_X_i(1)+t(1):DT:t(end)+time_conv_X_i(end)];

%Demodulating signal Y_Q_mod
Y_Q = Ts*conv(phi,Y_Q_demod);
time_of_X_q = [time_conv_X_q(1)+t(1):DT:t(end)+time_conv_X_q(end)];

% figure(9);
% subplot(2,1,1);
% plot(time_of_X_i,Y_I);
% xlabel('Time(s)');
% ylabel('Y_I');
% title('Convolution of demodulated Y_I and srrc pulse \phi');
% 
% subplot(2,1,2);
% plot(time_of_X_q,Y_Q);
% xlabel('Time(s)');
% ylabel('Y_Q');
% title('Convolution of demodulated Y_Q and srrc pulse \phi');

%Creating periodogram of Y__i
t_total_conv_X_i = time_of_X_i(end) - time_of_X_i(1);
P_Y_i_F = (abs(fftshift(fft(Y_I,Nf)*Ts)).^2)./(t_total_conv_X_i);

%Creating periodogram of Y__q
t_total_conv_X_q = time_of_X_q(end) - time_of_X_q(1);
P_Y_q_F = (abs(fftshift(fft(Y_Q,Nf)*Ts)).^2)./(t_total_conv_X_q);


% figure(10);
% subplot(2,1,1);
% plot(F,P_Y_i_F);
% xlabel('Frequency(Hz)');
% ylabel('P_{Y{_I}}_F');
% title('Periodogram of Y_I');
% 
% subplot(2,1,2);
% plot(F,P_Y_q_F);
% xlabel('Frequency(Hz)');
% ylabel('P_{Y{_Q}}_F');
% title('Periodogram of Y_Q');

%%%%%part 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_i = Y_I(2*A_pulse*over+1:over:2*A_pulse*over +N*over);
y_q =  Y_Q(2*A_pulse*over+1:over:2*A_pulse*over +N*over);

% figure(11);
% scatter(y_i,y_q);

%%%%%part 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimating symbol of Sequence x_i
est_X_i = detect_4_PAM(y_i,A);

est_X_q = detect_4_PAM(y_q,A);

%%%%%part 13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimating the propability of symbol error of Sequence X_I and X_Q

E_symbols_X_i = [est_X_i ~= X_I];

E_symbols_X_q = [est_X_q ~= X_Q];

E_symb = [E_symbols_X_i | E_symbols_X_q];
E_symbols = sum(E_symb);

%%%%%part 14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimating the propability of bit error of Sequence b_i and b_q

b_i = PAM_4_to_bits(est_X_i,A);

b_q = PAM_4_to_bits(est_X_q,A);

%%%%%part 15 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bits = [b_i ;b_q];

E_b = [b ~= bits];

E_bits = sum(E_b);
%------------------------PART B-----------------------%
%%%%%%%%%%%%%%%%%%%%%%Part 1%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_E_symbols = E_symbols/N;
P_E_bits = E_bits/(4*N);

Average_E_symbols = Average_E_symbols + P_E_symbols;
Average_E_bits = Average_E_bits + P_E_bits;

end

Tot_E_symbol(j) = Average_E_symbols/K;

Tot_E_bits(j) = Average_E_bits/K;

sigma_N = sqrt((Ts*(sigma^2))/2);

Theoretical_E_symbol = 3*Q((A/sigma_N)) - (9/4)*(Q((A/sigma_N)))^2;
Theoretical_E_bit = Theoretical_E_symbol/(log2(M));

Theoretical_Tot_E_symbol(j) = Theoretical_E_symbol;
Theoretical_Tot_E_bits(j) = Theoretical_E_bit;
end
%%%%%%%%%%%%%%%%%%%%%Part 2%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12);
semilogy(SNR_dB,Tot_E_symbol);
hold on;
semilogy(SNR_dB,Theoretical_Tot_E_symbol);
legend('Experimetnal SER','Theoretical SER');
xlabel('SNR in dB');
ylabel('Symbol Error Rate for 16-QAM');
title('Symbol Error Rate (SER) for SNR_{db} [0:2:16]');
hold off;

figure(13);
semilogy(SNR_dB,Tot_E_bits);
hold on;
semilogy(SNR_dB,Theoretical_Tot_E_bits);
legend('Experimetnal BER','Theoretical BER');
xlabel('SNR in dB');
ylabel('Bit Error Rate for 16-QAM');
title('Bit Error Rate (BER) for SNR_{db} [0:2:16]');
hold off;
