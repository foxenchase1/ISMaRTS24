% Calculus of Variations Brachristochrone

% Given functions for length of the brachristochrone and time to roll down
% the curve, we determine how the position of the endpoint (b,-1) where
% b\in[1,10] impacts the length and time.

close all
clear 
clc
format long

N = 50;
h = 10/N;
b1 = [1:h:10];
g = 9.81;
eps = 0.00001;
theta_0 = zeros(1,10);

time_solns = zeros(1,10);
length_solns = zeros(1,10);
line_len_solns = zeros(1,10);
line_time_solns = zeros(1,10);
parab_len_solns = zeros(1,10);
parab_time_solns = zeros(1,10);

parab_arc_fun = @(t,beta) sqrt(1+4*t.^2-4*t*beta+beta.^2);
parab_time_fun = @(t,beta) sqrt(1+(2*t-beta).^2)/sqrt(beta*t-t.^2);
arc_integrand_func = @(t) sqrt(1-cos(t));

syms theta

for i=1:length(b1)
    b = b1(i);
    bf = b == (theta-sin(theta))/(1-cos(theta));
    solx = vpasolve(bf,theta, [0 , 2*pi]);
    theta_0(i) = solx;
    a = b/(theta_0(i)-sin(theta_0(i)));
    beta1 = (1+b^2)/b;
    time_solns(i) = sqrt(b/(g*(theta_0(i)-sin(theta_0(i)))))*theta_0(i);
    %length_solns(i) = (-4*b/(theta_0(i)+sin(theta_0(i))))*(cos(theta_0(i)/2)-1);
    length_solns(i) = sqrt(2)*a*integral(arc_integrand_func,0,theta_0(i));
    line_len_solns(i) = sqrt(1+b^2);
    line_time_solns(i) = sqrt((2*b^2+2)/g);  
    parab_len_solns(i) = integral(@(t) parab_arc_fun(t,beta1), 0, b);
    parab_time_solns(i) = (1/sqrt(2*g))*integral(@(t) parab_time_fun(t,beta1), 0+eps, b-eps, 'ArrayValued',true);
end

figure
subplot(2,1,1)
plot(b1,length_solns,'r',b1,line_len_solns,'m',b1,parab_len_solns,'b')
xlabel('b Values')
title('Arc Lengths for Brach., Line, and Parab. to (b,-1)')
legend('Brach.','Linear','Quadratic','Location','northwest')
subplot(2,1,2)
plot(b1, time_solns,'r', b1, line_time_solns,'m', b1, parab_time_solns,'b')
xlabel('b Values')
title('Descent Times for Brach., Line, and Parab. to (b,-1)')
legend('Brach.','Linear','Quadratic','Location','northwest')
