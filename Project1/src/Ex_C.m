clear all;
close all;
T = 0.1;
over = 10;
a = 0.5;
A = 5;
N = 50;
Ts = T/over;
DT = Ts;
% part c
%part c1
b = (sign(randn(N,1)) + 1)/2;%The following of bits
%part c2
X = bits_to_2PAM(b);%The following of 2PAM symbols
%part b of c2
X_delta = (1/Ts)*upsample(X,over);%Upsampling of X
t = [0:Ts:(N*T - Ts)];
plot(t,X_delta,'LineWidth',1.5)
xlabel('Time of X_d');
ylabel('X_delta');
title('The signal X_d');
figure();
%part c of c2
%Modulation of the signal
[phi,t1] = srrc_pulse(T, Ts, A, a);
conv_phi_Xd = conv(phi,X_delta)*Ts;%Convolution of X_delta and phi
time_of_conv = [(t1(1)+t(1)):Ts:(t1(end) + t(end))];%Time of convolution of X_delta and phi
plot(time_of_conv,conv_phi_Xd);
xlabel('Time of convolution');
ylabel('Convolution of phi and X_d');
title('The convolution of phi and X_d');
figure();
%part d of c2
%Demodulation of the signal
inv_phi = phi;%The signal phi is even
Z = conv(phi,conv_phi_Xd)*Ts;%Z(t)
time_of_Z = [t1(1) + time_of_conv(1):Ts:t1(end)+time_of_conv(end)];%Time of Z(t)
plot(time_of_Z,Z)xlabel('Time of Z');
ylabel('Z(t)');
title('The signal of Z(t)');
figure();
plot(time_of_Z,Z)
hold on;
stem([0:1:N-1]*T,X)
xlabel('Time of Z');
ylabel('Z(t)');
title('The signal of Z(t)');