load const_N.mat
%% stress balance checks out

params = results_constN.params;
for t=1:10:100
    plot_vars_new(results_constN,t,params)
end
h = results_constN.hs(:,1);
u = results_constN.us(:,1);
xg = results_constN.xgs(1);
sigma_elem = params.sigma_elem;
sigma_h = params.fixed_N_grid;
sigma = params.sigma; 
x_grid_ice = sigma*xg*params.x0;
h_interp = interp1(sigma_elem,h.*params.h0,sigma,"linear","extrap");
N_interp = interp1(sigma_h,params.N_scaled./params.N0,sigma.*xg,"linear","extrap");
shear_stress = params.C.*N_interp(2:end)*params.N0.*(u(2:end).*params.u0./(u(2:end).*params.u0+params.As*(params.C*N_interp(2:end)*params.N0).^params.n)).^(1/params.n); 
b = -bed_schoof(x_grid_ice,params);
driving_stress = params.rho_i*params.g*h_interp(2:end).*gradient(h_interp(2:end)-b(2:end))./gradient(x_grid_ice(2:end));
shear_and_driving = shear_stress + driving_stress;
figure()
dudx = gradient(u(2:end).*params.u0)./gradient(x_grid_ice(2:end));
long_stress = gradient(2*params.A^(-1/params.n)*h_interp(2:end).*abs(dudx).^(1/params.n -1).*dudx)./gradient(x_grid_ice(2:end));
plot(x_grid_ice(2:end)./1000,shear_stress,'DisplayName', 'Shear stress');
hold on;
plot(x_grid_ice(2:end)./1000,driving_stress,'DisplayName','Driving stress');
plot(x_grid_ice(2:end)./1000,long_stress,'DisplayName', 'Longitudinal stress');
plot(x_grid_ice(2:end)./1000,shear_and_driving,'DisplayName', 'Shear + driving stress');
legend;
xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
title('$\frac{\partial }{\partial x}\left[ 2\bar{A}^{-1/n}h_g\left|\frac{\partial u}{\partial x}\right|^{1/n-1}\frac{\partial u}{\partial x}\right] - CN(\frac{u}{u+A_sC^nN^n})^{1/n} -\rho_igh_g\frac{\partial (h_g-b)}{\partial x}=0$','Interpreter','latex')


%% check the basal shear stress term as a function of u

umin = min(min(results_constN.us));
umax = max(max(results_constN.us));
u_range = linspace(umin,umax);
u_mesh = ones(length(u_range)).*u_range;

N_range = linspace(0,max(N_interp));
N_mesh = (ones(length(N_range)).*N_range)';

nondim_tau = (params.gamma.*N_mesh.*(u_mesh./(u_mesh+N_mesh.^params.n)).^(1./params.n));
figure()
surface(u_mesh,N_mesh,nondim_tau,'EdgeColor','none');c = colorbar;xlabel('u');ylabel(c,'tau');ylabel('N');

load const_powerlaw.mat;
figure();
params = results_constN.params;
lambda_range = linspace(min(params.lambda),max(params.lambda));
lambda_mesh = (ones(length(lambda_range)).*lambda_range)';
nondim_tau = lambda_mesh.*u_mesh;
surface(u_mesh,lambda_mesh,nondim_tau,'EdgeColor','none');c = colorbar;xlabel('u');ylabel(c,'tau');ylabel('lambda');


function plot_vars_new(results,time,params)
    sigma_elem = params.sigma_elem;
    sigma_h = params.sigma_h;
    sigma = params.sigma; 
    h = results.hs(:,time);
    u = results.us(:,time);
    xg = results.xgs(time);
        
   subplot(2,1,1);
    plot(sigma_elem.*xg.*params.x0./1000,h);ylabel('h');hold on;
    if time==81
        plot(linspace(0,1)*1500,bed_schoof(linspace(0,1).*1500e3,params)./params.h0,'-k'); 
    end
    xlim([1300 1375]);

    ax5 = subplot(2,1,2);
    plot(sigma.*xg.*params.x0./1000,u);hold on;ylabel('u');xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    xlim([1300 1375]);

    

    
end

