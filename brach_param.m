% Brachistochrone Parameter Function
function [x_brach, y_brach] = brach_param(x)
% Input vector x = [theta a]
% Output vector param = [x y]
x_brach = x(2)*(x(1)-sin(x(1)));
y_brach = x(2)*(1-cos(x(1)));
return [x_brach, y_brach]