%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Part B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part B4------------------------------------------------------------------------------
T = 0.001;
over = 10;
Ts = T/over;
A = 4;
a = 0.5;
Fs = 1/Ts;
Nf = 2048;
DT = 1/Fs;N = 100;
F = [-Fs/2:Fs/Nf:Fs/2 - Fs/Nf];
K = 1000;
F0 = 1/T;
[phi,t_phi] = srrc_pulse(T, Ts, A, a);
y1 = abs(fftshift(fft(phi,Nf)*DT)).^2;
Theta = 2*pi*randi(1,1);
sum_of_realizations = zeros(1,Nf);
for i = 1:K
b = (sign(randn(N,1)) + 1)/2;
X_n = bits_to_2PAM(b);
X_d = (1/Ts)*upsample(X_n,over);
t_Xd = [0:Ts:(N*T - Ts)];
X = conv(phi,X_d)*DT;
time_of_X = [(t_Xd(1)+t_phi(1)):DT:(t_Xd(end) + t_phi(end))];
Y = X.*cos(2*pi*F0.*time_of_X + Theta);
t_total = time_of_X(end)-time_of_X(1);
P_Y_F = (abs(fftshift(fft(Y,Nf)*Ts)).^2)./(t_total);
sum_of_realizations= sum_of_realizations +P_Y_F;
end
average_sum_of_realizations = sum_of_realizations/K;
theoretical_P_Y_F = (1/T).*y1;
figure(11);
semilogy(F,P_Y_F);
xlabel('Frequency(Hz)','FontSize',5);
ylabel('P_Y(F)','FontSize',5);
title('Periodogram an implementation of Y with 2PAM in logarithmic y axis ');
figure(12);semilogy(F,average_sum_of_realizations);
xlabel('Frequency(Hz)','FontSize',10);
ylabel('P_Y(F)','FontSize',10);
title('Power spectral density of Y');