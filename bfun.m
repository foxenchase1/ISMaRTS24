function cycl1 = bfun(theta, b)
% The zeros of this function w.r.t. theta is the corresponding 
% theta_0 for the parameterization of x=b
cycl1 = b*cos(theta)-sin(theta)+theta-b;