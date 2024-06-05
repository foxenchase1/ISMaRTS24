function isagb_time_integral = isg_t_integ(r,Rmax,R_apse)
% This function computes the descent time for the brachistochrone
% from a specified maximum radius to a minimum radius under inverse square
% attractive gravity (radius is to center of attracting body).
if Rmax <= R_apse
    disp('ERROR! You must have Rmax > R_apse.')
else
    prelim_func1 = @(zeta) (1/2)*(Rmax+R_apse+(Rmax-R_apse)*cos(zeta));
    integrand_isag_time = @(zeta) (prelim_func1(zeta))/sqrt(2*((prelim_func1(zeta))+R_apse*prelim_func1(zeta)+R_apse^2/(Rmax-R_apse)));
    upper_bd_func = @(rad) acos((2*rad-Rmax-R_apse)/(Rmax-R_apse));
    isagb_time_integral = zeros(1,length(r));
    upper_bds = zeros(1,length(r));
    for j=1:length(r)
        upper_bd = upper_bd_func(r(j));
        upper_bds(j) = upper_bd;
        isagb_time_integral_comp = integral(integrand_isag_time,0,upper_bd,'ArrayValued',true);
        isagb_time_integral(j) = isagb_time_integral_comp;
    end
end