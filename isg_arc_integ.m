function isagb_arc_integral = isg_arc_integ(r,Rmax,R_apse);
% This function computes the arc length for the brachistochrone
% from a specified maximum radius to a minimum radius under inverse square
% attractive gravity (radius is to center of attracting body).

if Rmax <= R_apse
    disp('ERROR! You must have Rmax > R_apse.')
else
    x = @(y) y^2+R_apse;
    arc_integrand = @(y) 2*x(y)/sqrt(x(y)+R_apse+(R_apse^2/(Rmax-R_apse)*x(y)));
    upper_bd = sqrt(Rmax-R_apse);
    lower_bd_func = @(r) sqrt(r-R_apse);
    lower_bds = zeros(1,length(r));
    isagb_arc_integral = zeros(1,length(r));
    for j=1:length(r)
        lower_bd = lower_bd_func(r(j));
        lower_bds(j) = lower_bd;
        isagb_arc_integral_comp = integral(arc_integrand,lower_bd,upper_bd,'ArrayValued',true);
        isagb_arc_integral(j) = isagb_arc_integral_comp;
    end
end