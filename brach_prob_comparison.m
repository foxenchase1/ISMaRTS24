%% Brachistochrone Solution Comparison
% Chase Foxen
% 4/1/24

% This program will plot a cycloid, parabola, and line from 0 to a
% specified x = b for every integer between 1 and b, inclusive.
% Two plots are then displayed which compare the descent times for each
% curve and the arc lengths of each curve.

clc
clear
close all 

b = input('What is the value of b? ');
b_check = isempty(b);
h = input('What step size? ');
h_check = isempty(h);
if b_check == 1
    disp('ERROR! You must specify a value of b.')
elseif h_check == 1
    disp('ERROR! You must specify a value of h.')
else
b1 = [1:1:b];
b2 = [1:h:b];
y = 1;
x = b;
x2 = [0:h:b];

g = 9.81;

func1 = @(x,t) (t-sin(t))/(1-cos(t))-x;
func2 = @(t) func1(x,t);
func3 = @(a,t) a*t-sin(t);
func4 = @(a,t) a-a*cos(t);

x_param = @(a,t) a*t-a*sin(t);
y_param = @(a,t) a-a*cos(t);
parab_func = @(b,x) x.^2-((1+b^2)/b)*x;
parab_arc_length_fun = @(b,x) sqrt(1+4*x.^2-4*x*((1+b^2)/b)+((1+b^2)/b)^2);
parab_time_fun = @(b,x) sqrt(1+(2.*x-((1+b^2)/b)).^2)/sqrt(((1+b^2)/b).*x-x.^2);
linear_func = @(b,x) (-1/b)*x;

t_last = fzero(func2,0.125);
t_fin_vector = [0:h:t_last];
a = x/(t_last-sin(t_last));

%descent_time = sqrt(b/(g*t-g*sin(t)))*t;
%arc_length = (-4*x*(cos(t/2)-1))/(t+sin(t));
%integrand_fun = @(t) sqrt(1-cos(t));
%arc_length_check = sqrt(2)*a*integral(integrand_fun,0,t_var);

% Memory for Plots and Other Data
parabs = zeros(b,length(x2));
lines = zeros(b,length(x2));

index_y_parab = zeros(1,b);
index_y_linear = zeros(1,b);

x3 = zeros(b,length(x2));
x_s = zeros(b,length(t_fin_vector));
y_s = zeros(b,length(t_fin_vector));
ts = zeros(1,b);
as = zeros(1,b);

parab_arc_lengths = zeros(1,b);
parab_descent_times = zeros(1,b);
lin_arc_lengths = zeros(1,b);
lin_descent_times = zeros(1,b);
cyc_arc_lengths = zeros(1,b);
cyc_descent_times = zeros(1,b);

for j=1:length(b1)
    b = b1(j);    
    x3(j,:) = x2;
    %Parabolic Path
    parab_1 = @(x) parab_func(b,x);
    parab_val = parab_1(x2);
    index_y_parab(j) = (b/h)+1;
    parabs(j,:) = parab_val;
    parabs(j,index_y_parab(j)+1:length(parab_val)) = 0;
    %Descent Times and Arc Lengths for Parabola
    parab_arc_integrand = @(x) parab_arc_length_fun(b,x);
    parab_arc_length = integral(parab_arc_integrand, 0+eps, b-eps);
    parab_arc_lengths(j) = parab_arc_length;
    parab_descent_integrand = @(x) parab_time_fun(b,x);
    parab_descent_time = (1/sqrt(2*g))*integral(parab_descent_integrand, 0+eps, b-eps,'ArrayValued',true);
    parab_descent_times(j) = parab_descent_time;
    %Linear Path
    line_1 = @(x) linear_func(b,x);
    linear_val = line_1(x2);
    index_y_linear(j) = (b/h)+1; 
    lines(j,:) = linear_val;
    lines(j,index_y_linear(j)+1:length(linear_val)) = 0;
    %Descent Times and Arc Lengths for Lines
    lin_arc_length = sqrt(1+b^2);
    lin_arc_lengths(j) = lin_arc_length;
    lin_descent_time = sqrt((2*b^2+2)/g);
    lin_descent_times(j) = lin_descent_time;
    %Cycloidal Path
    func2 = @(t) func1(b,t);
    t = fzero(func2,0.125);
    ts(j) = t;
    a = b/(t-sin(t));
    as(j) = a;
    t1 = [0:h:t];
    x1 = @(t) x_param(a,t);
    y1 = @(t) y_param(a,t);
    x_val = x1(t1);
    y_val = y1(t1);
    x_s(j,:) = [x_val,zeros(1,length(t_fin_vector)-length(x_val))];
    y_s(j,:) = [y_val,zeros(1,length(t_fin_vector)-length(y_val))];
    %Descent Times and Arc Lengths for Brach.
    cyc_descent_time = sqrt(b/(g*t-g*sin(t)))*t;
    cyc_descent_times(j) = cyc_descent_time;
    integrand_fun = @(t) sqrt(1-cos(t));
    cyc_arc_length = sqrt(2)*a*integral(integrand_fun,0,t);
    cyc_arc_lengths(j) = cyc_arc_length;
end

if b < 20
    figure(1)
    for j = 1:length(b1)
        plot(x3(j,1:index_y_parab(j)),parabs(j,1:index_y_parab(j)),x3(j,1:index_y_linear(j)),lines(j,1:index_y_linear(j)),x_s(j,:),-y_s(j,:)); hold on
    end
    xlabel('b Values')
    ylabel('y')
    title('Cycloidal, Linear, and Parabolic Paths to (b,-1)')

    figure(2)
    plot(b1,cyc_descent_times,'r',b1,parab_descent_times,'m',b1,lin_descent_times,'b')
    title('Descent Times of Curves')
    xlabel('b Value')
    ylabel('Descent Time of Curve')
    legend('Cycloid','Parabola','Line',Location='northwest')

    figure(3)
    plot(b1,cyc_arc_lengths,'r',b1,parab_arc_lengths,'m',b1,lin_arc_lengths,'b')
    title('Arc Lengths of Curves')
    xlabel('b Value')
    ylabel('Arc Length of Curve')
    legend('Cycloid','Parabola','Line',Location='northwest')
elseif b < 50
    figure(1)
    plot(b1,cyc_descent_times,'r',b1,parab_descent_times,'m',b1,lin_descent_times,'b')
    title('Descent Times of Curves')
    xlabel('b Value')
    ylabel('Descent Time of Curve')
    legend('Cycloid','Parabola','Line',Location='northwest')

    figure(2)
    plot(b1,cyc_arc_lengths,'r',b1,parab_arc_lengths,'m',b1,lin_arc_lengths,'b')
    title('Arc Lengths of Curves')
    xlabel('b Value')
    ylabel('Arc Length of Curve')
    legend('Cycloid','Parabola','Line',Location='northwest')
else
    figure(1)
    plot(b1,cyc_descent_times,'r',b1,parab_descent_times,'m',b1,lin_descent_times,'b')
    title('Descent Times of Curves')
    xlabel('b Value')
    ylabel('Descent Time of Curve')
    legend('Cycloid','Parabola','Line',Location='northwest')
end
end