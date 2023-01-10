close all;
clear all;
T = 0.001;
over = 10;
Ts = T/over;
A = 4;
a = 0.5;
Fs = 1/Ts;
Nf = 2048;
DT = 1/Fs;
N = 100;
F = [-Fs/2:Fs/Nf:Fs/2 - Fs/Nf];
K = 1000;
%Part A1-------------------------------------------------
[phi,t_phi] = srrc_pulse(T, Ts, A, a);%Creating srrc pulsess
plot(t_phi,phi);
figure(1);
y1 = abs(fftshift(fft(phi,Nf)*DT)).^2;%Creating energy spectral density
semilogy(F,y1);xlabel('Frequency(Hz)','FontSize',10);
ylabel('|\Phi(F)|^2','FontSize',10);
title('Energy Spectral Density of \Phi','FontSize',10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_of_realizations = zeros(1,Nf);
%Part A2-------------------------------------------------
for i = 1:K
b = (sign(randn(N,1)) + 1)/2;
X_n = bits_to_2PAM(b);
X_d = (1/Ts)*upsample(X_n,over);
t_Xd = [0:Ts:(N*T - Ts)];
X = conv(phi,X_d)*DT;
time_of_X = [(t_Xd(1)+t_phi(1)):DT:(t_Xd(end) + t_phi(end))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part A3
t_total = time_of_X(end)-time_of_X(1);
P_X_F = (abs(fftshift(fft(X,Nf)*Ts)).^2)./(t_total);
sum_of_realizations = sum_of_realizations +P_X_F;
end
average_sum_of_realizations = sum_of_realizations/K;
theoretical_P_X_F = (1/T).*y1;
figure(2);
plot(F,P_X_F);
xlabel('Frequency(Hz)','FontSize',10);
ylabel('P_X(F)','FontSize',10);
title('Periodogram of an implementation of X with bit 2PAM ');
figure(3);
semilogy(F,P_X_F);
xlabel('Frequency(Hz)','FontSize',5);
ylabel('P_X(F)','FontSize',5);
title('Periodogram of an implementation of X with bit 2PAM in logarithmic y axis');
figure(4);semilogy(F,average_sum_of_realizations);
hold on;
semilogy(F,theoretical_P_X_F);
xlabel('Frequency(Hz)','FontSize',10);
ylabel('P_X(F)','FontSize',10);
title('Power spectral density of X');
hold off;
%Part A4------------------------------------------------------------------------------
sum_of_realizations = zeros(1,Nf);
for i = 1:K
b = (sign(randn(N,1)) + 1)/2;
X_n = bits_to_4PAM(b);
X_d = (1/Ts)*upsample(X_n,over);
t_Xd = [0:Ts:((N*T)/2 - Ts)];
X = conv(phi,X_d)*DT;
time_of_X = [(t_Xd(1)+t_phi(1)):DT:(t_Xd(end) + t_phi(end))];
t_total = time_of_X(end)-time_of_X(1);
P_X_F = (abs(fftshift(fft(X,Nf)*Ts)).^2)./t_total;
sum_of_realizations = sum_of_realizations +P_X_F;
end
average_sum_of_realizations = sum_of_realizations/K;
theoretical_P_X_F = (5/T).*y1;
figure(5);
plot(F,P_X_F);
xlabel('Frequency(Hz)','FontSize',10);
ylabel('P_X(F)','FontSize',10);
title('Periodogram of an implementation of X with bit 4PAM ');
figure(6);
semilogy(F,P_X_F);
xlabel('Frequency(Hz)','FontSize',5);
ylabel('P_X(F)','FontSize',5);title('Periodogram of an implementation of X with bit 4PAM in logarithmic y axis');
figure(7);
semilogy(F,average_sum_of_realizations);
hold on;
semilogy(F,theoretical_P_X_F);
xlabel('Frequency(Hz)','FontSize',10);
ylabel('P_X(F)','FontSize',10);
title('Power spectral density of X');
hold off;
%Part A5------------------------------------------------------------------------------
new_over = 2*over;
new_T = 2*T;
Nf = 4096;
new_F = [-Fs/2:Fs/Nf:Fs/2 - Fs/Nf];
[phi,t_phi] = srrc_pulse(new_T, Ts, A, a);
y1 = abs(fftshift(fft(phi,Nf)*DT)).^2;
sum_of_realizations = zeros(1,Nf);
P_X_F = zeros(1,Nf);
average_sum_of_realizations = zeros(1,Nf);
for i = 1:K
b = (sign(randn(N,1)) + 1)/2;
X_n = bits_to_2PAM(b);
X_d = (1/Ts)*upsample(X_n,new_over);
t_Xd = [0:Ts:(N*new_T - Ts)];
X = conv(phi,X_d)*DT;
time_of_X = [(t_Xd(1)+t_phi(1)):DT:(t_Xd(end) + t_phi(end))];
t_total = time_of_X(end)-time_of_X(1);
P_X_F = (abs(fftshift(fft(X,Nf)*Ts)).^2)./t_total;
sum_of_realizations = sum_of_realizations +P_X_F;end
average_sum_of_realizations = sum_of_realizations/K;
theoretical_P_X_F = (1/new_T).*y1;
figure(8);
plot(new_F,P_X_F);
xlabel('Frequency(Hz)','FontSize',5);
ylabel('P_X(F)','FontSize',5);
title('Periodogram of an implementation of X with bit 2PAM ,new_T = 2T');
figure(9);
semilogy(new_F,P_X_F);
xlabel('Frequency(Hz)','FontSize',5);
ylabel('P_X(F)','FontSize',5);
title('Periodogram of an implementation of X with bit 2PAM in logarithmic y axis, new_T =2T');
figure(10);
semilogy(new_F,average_sum_of_realizations);
hold on;
semilogy(new_F,theoretical_P_X_F);
xlabel('Frequency(Hz)','FontSize',10);
ylabel('P_X(F)','FontSize',10);
title('Power spectral density of X, new_T = 2T');
hold off;