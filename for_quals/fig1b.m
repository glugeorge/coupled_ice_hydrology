load quad_oneway.mat;
figure();
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 
% For oneway
    xg = params.xg;
    h = params.h;
    u = params.u;
    Q = Qs(end,:)';
    N = Ns(end,:)';
    S = Ss(end,:)';

% end
x_grid_hydro = sigma_h*xg*params.x0;
h_interp = interp1(sigma_elem,h.*params.h0,sigma_h,"linear","extrap");
u_interp = interp1(sigma,u.*params.u0,sigma_h,"linear","extrap");
b = -bed(x_grid_hydro,params);
psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);

subplot(5,1,1);
plot(x_grid_hydro./1000,params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro),'DisplayName','$\rho_wg\frac{\partial b}{\partial x}$');
hold on; 
plot(x_grid_hydro./1000,gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro),'DisplayName','$\rho_ig\frac{\partial h}{\partial x}$')
plot(x_grid_hydro./1000,psi,'DisplayName', '$\psi$');
legend('Interpreter','latex');
title('$\psi = \rho_wg\frac{\partial b}{\partial x} - \rho_ig\frac{\partial h}{\partial x}$','Interpreter','latex');

m_over_rho_i = params.m0.*abs(Q).^3./(params.rho_i.*S.^(8/3)); % from m'=|Q'|^3/S'^8/3
KSN3 = params.K0.*S.*params.S0.*(N.*params.N0).^3;
u_advect = u_interp.*gradient(S.*params.S0)./gradient(x_grid_hydro);

subplot(5,1,2);
plot(x_grid_hydro./1000,m_over_rho_i,'k','DisplayName', '$\frac{m}{\rho_i}$');
hold on;
plot(x_grid_hydro./1000,KSN3,'--','DisplayName', '$K_0SN^3$');
plot(x_grid_hydro./1000,u_advect,':','DisplayName', '$u\frac{\partial S}{\partial x}$');
legend('Interpreter','latex');
title('$\frac{\partial S}{\partial t_h} = \frac{m}{\rho_i}-K_0SN^3 - u\frac{\partial S}{\partial x}$',Interpreter='latex')

subplot(5,1,3)
dQdx = gradient(Q.*params.Q0)./gradient(x_grid_hydro);
m_over_rho_w = m_over_rho_i.*params.rho_i./params.rho_w;
M = ones(size(x_grid_hydro)).*params.M*params.M0;
plot(x_grid_hydro./1000,dQdx,'DisplayName', '$\frac{\partial Q}{\partial x}$');
hold on;
plot(x_grid_hydro./1000,m_over_rho_w,'DisplayName', '$\frac{m}{\rho_w}$');
plot(x_grid_hydro./1000,M,'DisplayName', '$M$');
legend('Interpreter','latex');
title('$\frac{\partial S}{\partial t_h} + \frac{\partial Q}{\partial x} = \frac{m}{\rho_w}+M$','Interpreter','latex')

subplot(5,1,4);
dNdx = gradient(N.*params.N0)./gradient(x_grid_hydro);
f_Q_S = params.f*params.rho_w*params.g.*Q.*params.Q0.*abs(Q.*params.Q0)./(S.*params.S0).^(8/3);
h_interp = interp1(sigma_elem,h,sigma_h,'linear','extrap');
b = -bed(x_grid_hydro,params);
psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);
plot(x_grid_hydro./1000,psi,'DisplayName', '$\psi$');
hold on;
plot(x_grid_hydro./1000,dNdx,'DisplayName', '$\frac{\partial N}{\partial x}$');
plot(x_grid_hydro./1000,f_Q_S,'DisplayName', '$f\rho_wg\frac{Q|Q|}{S^{8/3}}$');

legend('Interpreter','latex');
title('$\psi+ \frac{\partial N}{\partial x} = f\rho_wg\frac{Q|Q|}{S^{8/3}}$','Interpreter','latex')

subplot(5,1,5);
x_grid_ice = sigma*xg*params.x0;
h_interp = interp1(sigma_elem,h.*params.h0,sigma,"linear","extrap");
N_interp = interp1(sigma_h,N,sigma,"linear","extrap");
shear_stress = params.C.*N_interp(2:end)*params.N0.*(u(2:end).*params.u0./(u(2:end).*params.u0+params.As*(params.C*N_interp(2:end)*params.N0).^params.n)).^(1/params.n); 
b = -bed(x_grid_ice,params);
driving_stress = params.rho_i*params.g*h_interp(2:end).*gradient(h_interp(2:end)-b(2:end))./gradient(x_grid_ice(2:end));
shear_and_driving = shear_stress + driving_stress;

dudx = gradient(u(2:end).*params.u0)./gradient(x_grid_ice(2:end));
long_stress = gradient(2*params.A^(-1/params.n)*h_interp(2:end).*abs(dudx).^(1/params.n -1).*dudx)./gradient(x_grid_ice(2:end));
plot(x_grid_ice(2:end)./1000,shear_stress,'DisplayName', 'Shear stress');
hold on;
plot(x_grid_ice(2:end)./1000,driving_stress,'DisplayName', 'Driving stress');
plot(x_grid_ice(2:end)./1000,long_stress,'DisplayName', 'Longitudinal stress');
legend;
xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
title('$\frac{\partial }{\partial x}\left[ 2\bar{A}^{-1/n}h_g\left|\frac{\partial u}{\partial x}\right|^{1/n-1}\frac{\partial u}{\partial x}\right] - CN(\frac{u}{u+A_sC^nN^n})^{1/n} -\rho_igh_g\frac{\partial (h_g-b)}{\partial x}=0$','Interpreter','latex')

    
%% Plot terms for coupled


function plot_terms(Q,N,S,h,u,xg,params)
    sigma_elem = params.sigma_elem;
    sigma_h = params.sigma_h;
    sigma = params.sigma; 
    
    % do ice plots on sigma_elem grid
    x_grid_ice = sigma*xg*params.x0;
    diff_ice = diff(x_grid_ice);
    figure();
    %% Plot 1 - ice mass conservation
    % Remember that since these are all in steady state, d/dt goes to zero
    accum = params.accum; % normalization done in flowline equations
    h_interp = interp1(sigma_elem,h.*params.h0,sigma,"linear","extrap");
    dhu_dx = gradient((h_interp.*u).*params.u0)./diff_ice(1);
    subplot(3,2,1);
    plot(x_grid_ice./1000,accum.*ones(size(x_grid_ice)),'DisplayName', 'Accumulation');
    hold on;
    plot(x_grid_ice./1000,dhu_dx,'--','DisplayName', '$\frac{\partial (hu)}{\partial x}$');
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\frac{\partial h_g}{\partial t} + \frac{\partial (h_gu)}{\partial x}= a$','Interpreter','latex')
    %% Plot 2: ice momentum conservation
    N_interp = interp1(sigma_h,N,sigma,"linear","extrap");
    shear_stress = params.shear_scale.*params.C.*N_interp(2:end)*params.N0.*(u(2:end).*params.u0./(u(2:end).*params.u0+params.As*(params.C*N_interp(2:end)*params.N0).^params.n)).^(1/params.n); 
    b = -bed(x_grid_ice,params);
    driving_stress = params.rho_i*params.g*h_interp(2:end).*gradient(h_interp(2:end)-b(2:end))./gradient(x_grid_ice(2:end));
    shear_and_driving = shear_stress + driving_stress;
    
    dudx = gradient(u(2:end).*params.u0)./gradient(x_grid_ice(2:end));
    long_stress = gradient(2*params.A^(-1/params.n)*h_interp(2:end).*abs(dudx).^(1/params.n -1).*dudx)./gradient(x_grid_ice(2:end));
    subplot(3,2,2);
    plot(x_grid_ice(2:end)./1000,shear_stress,'DisplayName', 'Shear stress');
    hold on;
    plot(x_grid_ice(2:end)./1000,driving_stress,'DisplayName', 'Driving stress');
    plot(x_grid_ice(2:end)./1000,long_stress,'DisplayName', 'Longitudinal stress');
    plot(x_grid_ice(2:end)./1000,shear_and_driving,'DisplayName', 'Shear + driving stress');
    legend;
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\frac{\partial }{\partial x}\left[ 2\bar{A}^{-1/n}h_g\left|\frac{\partial u}{\partial x}\right|^{1/n-1}\frac{\partial u}{\partial x}\right] - CN(\frac{u}{u+A_sC^nN^n})^{1/n} -\rho_igh_g\frac{\partial (h_g-b)}{\partial x}=0$','Interpreter','latex')

    
    %% Plot 3: Channel Cross sectional area
    x_grid_hydro = sigma_h*xg*params.x0;
    m_over_rho_i = params.m0.*abs(Q).^3./(params.rho_i.*S.^(8/3)); % from m'=|Q'|^3/S'^8/3
    KSN3 = params.K0.*S.*params.S0.*(N.*params.N0).^3;
    u_interp = interp1(sigma,u.*params.u0,sigma_h,"linear","extrap");
    u_advect = u_interp.*gradient(S.*params.S0)./gradient(x_grid_hydro);
    
    subplot(3,2,3);
    plot(x_grid_hydro./1000,m_over_rho_i,'k','DisplayName', 'Melt term');
    hold on;
    plot(x_grid_hydro./1000,KSN3,'--','DisplayName', '$K_0SN^3$');
    plot(x_grid_hydro./1000,u_advect,':','DisplayName', 'Advection term');
    plot(x_grid_hydro./1000,KSN3 + u_advect,'x','DisplayName', '$K_0SN^3$ + Advection term');
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\frac{\partial S}{\partial t_h} = \frac{m}{\rho_i}-K_0SN^3 - u\frac{\partial S}{\partial x}$',Interpreter='latex')
    %% Plot 4: Conservation of mass
    dQdx = gradient(Q.*params.Q0)./gradient(x_grid_hydro);
    m_over_rho_w = m_over_rho_i.*params.rho_i./params.rho_w;
    M = ones(size(x_grid_hydro)).*params.M*params.M0;
    subplot(3,2,4);
    plot(x_grid_hydro./1000,dQdx,'DisplayName', '$\frac{\partial Q}{\partial x}$');
    hold on;
    plot(x_grid_hydro./1000,m_over_rho_w,'DisplayName', 'Melt term');
    plot(x_grid_hydro./1000,M,'DisplayName', 'Supply term');
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\frac{\partial S}{\partial t_h} + \frac{\partial Q}{\partial x} = \frac{m}{\rho_w}+M$','Interpreter','latex')
    %% Plot 5: Momentum conservation
    dNdx = gradient(N.*params.N0)./gradient(x_grid_hydro);
    f_Q_S = params.f*params.rho_w*params.g.*Q.*params.Q0.*abs(Q.*params.Q0)./(S.*params.S0).^(8/3);
    h_interp = interp1(sigma_elem,h,sigma_h,'linear','extrap');
    b = -bed(x_grid_hydro,params);
    psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);
    subplot(3,2,5);
    plot(x_grid_hydro./1000,psi,'DisplayName', '$\psi$');
    hold on;
    plot(x_grid_hydro./1000,dNdx,'DisplayName', '$\frac{\partial N}{\partial x}$');
    plot(x_grid_hydro./1000,f_Q_S,'DisplayName', '$f\rho_wg\frac{Q|Q|}{S^{8/3}}$');
    plot(x_grid_hydro./1000,psi+dNdx,'--','DisplayName', '$\psi+\frac{\partial N}{\partial x}$');
    
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\psi+ \frac{\partial N}{\partial x} = f\rho_wg\frac{Q|Q|}{S^{8/3}}$','Interpreter','latex')
    %% Plot 6: Energy
    melt_term = m_over_rho_i.*params.rho_i*params.L;
    flow_through_potential = Q.*params.Q0 .*(psi + dNdx);
    subplot(3,2,6);
    plot(x_grid_hydro./1000,melt_term,'DisplayName', 'Melt');
    hold on;
    plot(x_grid_hydro./1000,flow_through_potential,'DisplayName', 'Flow through potential');
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$mL = Q(\psi + \frac{\partial N}{\partial x})$','Interpreter','latex')

    figure()
    plot(x_grid_hydro./1000,params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro),'DisplayName','$\rho_wg\frac{\partial b}{\partial x}$');
    hold on; 
    plot(x_grid_hydro./1000,gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro),'DisplayName','$\rho_ig\frac{\partial h}{\partial x}$')
    plot(x_grid_hydro./1000,psi,'DisplayName', '$\psi$');    
    legend('Interpreter','latex');


end 

