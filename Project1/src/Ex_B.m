clear all;
close all;
%------Second part of the project B-----%
T = 0.01;
over = 10;
Ts = T/over;
DT = Ts;
A = 5;
Fs = 1/Ts;
%----Part B----%
%---for roll-off a = 0---%
a = [0 0.5 1];
for i = 1 : length(a)
[phi,t] = srrc_pulse(T, Ts, A, a(i));
figure();
plot(t,phi);
xlabel('Time(s)');
ylabel('\phi(t)');
title(['The SRRC pulse for a = ',num2str(a(i))]);
axis([-0.05 0.05 -5 15]);
%%%%Part B1---------------------------------------------------
figure();
%k = 0
phi_initial_0 = phi;%change ph0 to phi1_delay_0
new_t_0 = t;
plot(new_t_0, phi_initial_0);
hold on;
% k = 1
%figure();
phi_initial_1 = [phi zeros(1,(T/Ts))];
phi_delay_1 = [zeros(1,(T/Ts)) phi];
new_t1 = [t [t(end)+ Ts:DT:t(end)+T]];
plot(new_t1,phi_delay_1);
hold on;
% k = 2
%figure();
phi_initial_2 = [phi zeros(1,2*(T/Ts))];
phi_delay_2 = [zeros(1,2*(T/Ts)) phi];
new_t2 = [t [t(end)+ Ts:DT:t(end)+2*T]];
plot(new_t2,phi_delay_2);hold on;
% k = 4
%figure();
phi_initial_4 = [phi zeros(1,4*(T/Ts))];
phi_delay_4 = [zeros(1,4*(T/Ts)) phi];
new_t4 = [t [t(end)+ Ts:DT:t(end)+4*T]];%I do not know if Ts is neseccesary
plot(new_t4,phi_delay_4);
hold on;
axis([new_t4(1) new_t4(end) min(phi) max(phi)]);
xlabel('Time(s)');
ylabel('SRRC pulses');
title(['Phi and its delayed versions of k for a =',num2str(a(i))]);
%%%%%%%Part B2------------------------------------------------
%mulitplication of the two signals for k = 0
phi_product_0 = (phi_initial_0.*phi_initial_0);%The product of the srrc pulse and its
delayed version for 0 T
figure();
plot(new_t_0,phi_product_0);
axis([new_t_0(1) new_t_0(end) min(phi_product_0) max(phi_product_0)]);
xlabel('Time(s)');
ylabel('phi\_product\_0');
title('The product of the srrc pulse and its delayed version for 0 T');
%mulitplication of the two signals for k =1
phi_product_1 = (phi_initial_1.*phi_delay_1);%The product of the srrc pulse and its delayed
version for 1 T
figure();
plot(new_t1,phi_product_1);
axis([new_t1(1) new_t1(end) min(phi_product_1) max(phi_product_1)]);
xlabel('Time(s)');
ylabel('phi\_product\_1');
title('The product of the srrc pulse and its delayed version for 1 T');
%mulitplication of the two signals for k =2
phi_product_2 = (phi_initial_2.*phi_delay_2);%The product of the srrc pulse and its delayed
version for 2 T
figure();
plot(new_t2,phi_product_2);
axis([new_t2(1) new_t2(end) min(phi_product_2) max(phi_product_2)]);
xlabel('Time(s)');
ylabel('phi\_product\_2');
title('The product of the srrc pulse and its delayed version for 2 T');
%mulitplication of the two signals for k = 4
phi_product_4 = (phi_initial_4.*phi_delay_4);%The product of the srrc pulse and its delayed
version for 4 T
figure();
plot(new_t4,phi_product_4);
axis([new_t4(1) new_t4(end) min(phi_product_4) max(phi_product_4)]);
xlabel('Time(s)');
ylabel('phi\_product\_4');
title('The product of the srrc pulse and its delayed version for 4 T');
%part B3------------------------------------------------------
%Calculate the integer of the product of srrc pulse and its delayed version
phi_sum_of_products_0 = sum(phi_product_0*Ts)phi_sum_of_products_1 = sum(phi_product_1*Ts)
phi_sum_of_products_2 = sum(phi_product_2*Ts)
phi_sum_of_products_4 = sum(phi_product_4*Ts)
end