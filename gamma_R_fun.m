function  gamma_R = gamma_R_fun(Rmax, Rmin, Theta_min, R_apse);
% An intermediate function for determining the apsidial distance R of 
% a bradchistochrone through two points in polar coordinates. 
%
% This is Equation (88) from the paper Isochrones and Brachistochrones
% by Garry J. Tee from the University of Auckland.

L_func = @(R_apse) 1/2*acos(((2*R_apse/Rmin)-Rmax-R_apse)/(Rmax-R_apse));
zed_func2 = @(lambda) (cos(lambda)).^2+Rmin*(sin(lambda)).^2;
gamma_R_integrand = @(lambda) 2*(zed_func2(lambda)-R_apse)/sqrt(zed_func2(lambda).^2+(Rmax-R_apse)*zed_func2(lambda)+Rmax-R_apse);
lower_bd = L_func(R_apse);
Aps_Theta = isg_ang_integ(R_apse,Rmax,R_apse);

if Aps_Theta >= Theta_min
    gamma_R = integral(gamma_R_integrand,lower_bd,pi/2,'ArrayValued',true);
else
    gamma_R = 2*Aps_Theta-integral(gamma_R_integrand,lower_bd,pi/2,'ArrayValued',true);
end

