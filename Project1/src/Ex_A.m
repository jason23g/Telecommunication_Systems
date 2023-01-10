clear all;
close all;
%------First part of the project A-----%
T = 0.01;
over = 10;
Ts = T/over;
A = 4;
%----Part A1----%
%Create the srrc pulses for different roll off factors
%---for roll-off a = 0---%
a1 = 0;
[phi1,t1] = srrc_pulse(T, Ts, A, a1);
plot(t1,phi1);
hold on;
%---for roll-off a = 0.5---%
a2 = 0.5;
[phi2,t2] = srrc_pulse(T, Ts, A, a2);
plot(t2,phi2);
hold on;
%---for roll-off a = 1---%
a3 = 1;
[phi3,t3] = srrc_pulse(T, Ts, A, a3);
grid on;
plot(t3,phi3);
xlabel('Time(s)','FontSize',10);
ylabel('SRRC(t)','FontSize',10);
title('SRRC pulses for differents roll off factors','FontSize',15);
axis([-0.04 0.04 -3 15]);
%----Part A2----%
figure();
%Calculation of the srrc pulses' fourier transform
Fs = 1/Ts;
Nf = 1024;
DT = 1/Fs;
F = [-Fs/2:Fs/Nf:Fs/2 - Fs/Nf];%The frequency spectrum
y1 = abs(fftshift(fft(phi1,Nf)*DT)).^2;%Create the energy spectral density for phi1
plot(F,y1);
hold on;
y2 = abs(fftshift(fft(phi2,Nf)*DT)).^2;%Create the energy spectral density for phi2
plot(F,y2);
hold on;
y3 = abs(fftshift(fft(phi3,Nf)*DT)).^2;%Create the energy spectral density for phi3
plot(F,y3);
grid on;
xlabel('Frequency(Hz)','FontSize',10);
ylabel('|\Phi(F)|^2','FontSize',10);
title('Energy Spectral Density','FontSize',15);
figure();
semilogy(F,y1);
hold on;
semilogy(F,y2);
hold on;
semilogy(F,y3);
xlabel('Frequency(Hz)','FontSize',10);
ylabel('|\Phi(F)|^2','FontSize',10);
title('Energy Spectral Density(logarithmic y axis)','FontSize',15);
%----Part A3----%BW1 = (1+a1)/(2*T) %spectrum of phi1
BW2 = (1+a2)/(2*T) %spectrum of phi2
BW3 = (1+a3)/(2*T) %spectrum of phi3
colorstring = 'kb';%k is for black color line and b for blue color line
c = T/(10^3);
plot(F,c,'LineWidth',2,'color',colorstring(1));%This line has black color
hold on;
c1 = T/(10^5);
plot(F,c1,'LineWidth',2,'color',colorstring(2));%This line has blue color