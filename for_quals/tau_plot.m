load coupled_init_cond.mat
u_interp = interp1(params.sigma,u.*params.u0,params.sigma_elem,"linear","extrap");
u_interp(1) = 0;
N = interp1(params.sigma_h,N,params.sigma_elem,"linear","extrap");
shear_stress = params.C.*N*params.N0.*(u_interp./(u_interp+params.As*(params.C*N*params.N0).^params.n)).^(1/params.n); 
b = -bed(xg.*params.sigma_elem.*params.x0,params);
driving_stress = params.rho_i*params.g*h*params.h0.*gradient(h*params.h0-b)./gradient(params.sigma_elem*xg*params.x0);

shear_and_driving = shear_stress + driving_stress;

dudx = gradient(u_interp)./gradient(params.sigma_elem*xg*params.x0);
long_stress = gradient(2*params.A^(-1/params.n)*h*params.h0.*abs(dudx).^(1/params.n -1).*dudx)./gradient(params.sigma_elem*xg*params.x0);

plot(params.sigma_elem,shear_stress,'DisplayName','Shear stress');
hold on;
plot(params.sigma_elem,driving_stress,'DisplayName','Driving stress');
plot(params.sigma_elem,long_stress,'DisplayName','Longitudinal stress');
plot(params.sigma_elem,shear_and_driving,'DisplayName','Shear + driving stress');
legend;

function b = bed(x,params)
    xsill = x>params.sill_min & x<params.sill_max;
    xdsill = x>=params.sill_max;
    sill_length = params.sill_max-params.sill_min;
    
    b = params.b0 + params.bx.*x;
    
    b(xsill) = params.b0 + (params.bx*params.sill_min) + params.sill_slope.*(x(xsill)-params.sill_min);

    b(xdsill) = params.b0 + (params.bx*params.sill_min) + params.sill_slope.*sill_length + ...
            params.bx*(x(xdsill)-params.sill_max);

end