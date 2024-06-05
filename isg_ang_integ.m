function isagb_angle_integral = isg_ang_integ(r,Rmax,R_apse)
% This function computes the angle for the brachistochrone
% from a specified maximum radius to a minimum radius under inverse square
% attractive gravity (radius is to center of attracting body).
if Rmax <= R_apse
    disp('ERROR! You must have Rmax > R_apse.')
else
    zed_func1 = @(rad) R_apse/rad;
    zed_func2 = @(lambda) (cos(lambda))^2+R_apse*(sin(lambda))^2;
    lambda_func = @(zed) 1/2*acos((2*zed-Rmax-R_apse)/(Rmax-R_apse));
    isg_ang_integrand = @(lambda) 2*(zed_func2(lambda)-R_apse)/sqrt((zed_func2(lambda))^2+(Rmax-R_apse)*zed_func2(lambda)+Rmax-R_apse);
    isagb_angle_integral = zeros(1,length(r));
    lower_bds = zeros(1,length(r));
    for j = 1:length(r)
        lower_bd = lambda_func(zed_func1(r(j)));
        lower_bds(j) = lower_bd;
        isagb_angle_integral_comp = integral(isg_ang_integrand,lower_bd,pi/2,'ArrayValued',true);
        isagb_angle_integral(j) = isagb_angle_integral_comp;
    end
end