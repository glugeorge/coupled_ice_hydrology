close all; clear;
load C_0.1_A_4.9_c.mat;
QNShuxg = results.init_cond;
params = results.params;
plot_terms(QNShuxg,params,'Initial')
QNShuxg = results.steady_state;
params = results.params;
plot_terms(QNShuxg,params,'Final')
function plot_terms(QNShuxg,params,I_F)
    sigma_elem = params.sigma_elem;
    sigma_h = params.sigma_h;
    sigma = params.sigma; 
    Q = QNShuxg(1:params.Nx);
    N = QNShuxg(params.Nx+1:2*params.Nx);
    S = QNShuxg(2*params.Nx+1:3*params.Nx);
    h = QNShuxg(3*params.Nx+1:4*params.Nx);
    u = QNShuxg(4*params.Nx+1:5*params.Nx);
    xg = QNShuxg(5*params.Nx+1);
    
    % do ice plots on sigma_elem grid
    x_grid_ice = sigma_elem*xg*params.x0;
    figure();
    %% Plot 1 - ice mass conservation
    % Remember that since these are all in steady state, d/dt goes to zero
    accum = params.accum; % normalization done in flowline equations
    u_interp = interp1(sigma,u.*params.u0,sigma_elem,"linear","extrap");
    dhu_dx = gradient(u_interp.*h.*params.h0)./gradient(x_grid_ice);
    subplot(3,2,1);
    plot(x_grid_ice./1000,accum.*ones(size(x_grid_ice)),'DisplayName', strcat(I_F,': ','Accumulation'));
    hold on;
    plot(x_grid_ice./1000,dhu_dx,'--','DisplayName', strcat(I_F,': ','$\frac{\partial (hu)}{\partial x}$'));
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\frac{\partial h_g}{\partial t} + \frac{\partial (h_gu)}{\partial x}= a$','Interpreter','latex')
    %% Plot 2: ice momentum conservation
    N_interp = interp1(sigma_h,N,sigma_elem,"linear","extrap");
    shear_stress = params.C.*N_interp(2:end)*params.N0.*(u_interp(2:end)./(u_interp(2:end)+params.As*(params.C*N_interp(2:end)*params.N0).^params.n)).^(1/params.n); 
    b = -bed(xg.*sigma_elem.*params.x0,params);
    driving_stress = params.rho_i*params.g*h(2:end)*params.h0.*gradient(h(2:end)*params.h0-b(2:end))./gradient(x_grid_ice(2:end));
    shear_and_driving = shear_stress + driving_stress;
    
    dudx = gradient(u_interp(2:end))./gradient(sigma_elem(2:end)*xg*params.x0);
    long_stress = gradient(2*params.A^(-1/params.n)*h(2:end)*params.h0.*abs(dudx).^(1/params.n -1).*dudx)./gradient(x_grid_ice(2:end));
    subplot(3,2,2);
    plot(x_grid_ice(2:end)./1000,shear_stress,'DisplayName', strcat(I_F,': ','Shear stress'));
    hold on;
    plot(x_grid_ice(2:end)./1000,driving_stress,'DisplayName', strcat(I_F,': ','Driving stress'));
    plot(x_grid_ice(2:end)./1000,long_stress,'DisplayName', strcat(I_F,': ','Longitudinal stress'));
    plot(x_grid_ice(2:end)./1000,shear_and_driving,'DisplayName', strcat(I_F,': ','Shear + driving stress'));
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
    plot(x_grid_hydro./1000,m_over_rho_i,'k','DisplayName', strcat(I_F,': ','Melt term'));
    hold on;
    plot(x_grid_hydro./1000,KSN3,'--','DisplayName', strcat(I_F,': ','$K_0SN^3$'));
    plot(x_grid_hydro./1000,u_advect,':','DisplayName', strcat(I_F,': ','Advection term'));
    plot(x_grid_hydro./1000,KSN3 + u_advect,'x','DisplayName', strcat(I_F,': ','$K_0SN^3$ + Advection term'));
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\frac{\partial S}{\partial t_h} = \frac{m}{\rho_i}-K_0SN^3 - u\frac{\partial S}{\partial x}$',Interpreter='latex')
    %% Plot 4: Conservation of mass
    dQdx = gradient(Q.*params.Q0)./gradient(x_grid_hydro);
    m_over_rho_w = m_over_rho_i.*params.rho_i./params.rho_w;
    M = ones(size(x_grid_hydro)).*params.M*params.M0;
    subplot(3,2,4);
    plot(x_grid_hydro./1000,dQdx,'DisplayName', strcat(I_F,': ','$\frac{\partial Q}{\partial x}$'));
    hold on;
    plot(x_grid_hydro./1000,m_over_rho_w,'DisplayName', strcat(I_F,': ','Melt term'));
    plot(x_grid_hydro./1000,M,'DisplayName', strcat(I_F,': ','Supply term'));
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\frac{\partial S}{\partial t_h} + \frac{\partial Q}{\partial x} = \frac{m}{\rho_w}+M$','Interpreter','latex')
    %% Plot 5: Momentum conservation
    dNdx = gradient(N.*params.N0)./gradient(x_grid_hydro);
    f_Q_S = params.f*params.rho_w*params.g.*Q.*params.Q0.*abs(Q.*params.Q0)./(S.*params.S0).^(8/3);
    h_interp = interp1(sigma_elem,h,sigma_h,'linear','extrap');
    psi = params.rho_w*params.g*sin(0.001)-gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);
    subplot(3,2,5);
    plot(x_grid_hydro./1000,psi,'DisplayName', strcat(I_F,': ','$\psi$'));
    hold on;
    plot(x_grid_hydro./1000,dNdx,'DisplayName', strcat(I_F,': ','$\frac{\partial N}{\partial x}$'));
    plot(x_grid_hydro./1000,f_Q_S,'DisplayName', strcat(I_F,': ','$f\rho_wg\frac{Q|Q|}{S^{8/3}}$'));
    plot(x_grid_hydro./1000,psi+dNdx,'--','DisplayName', strcat(I_F,': ','$\psi+\frac{\partial N}{\partial x}$'));
    
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\psi+ \frac{\partial N}{\partial x} = f\rho_wg\frac{Q|Q|}{S^{8/3}}$','Interpreter','latex')
    %% Plot 6: Energy
    melt_term = m_over_rho_i.*params.rho_i*params.L;
    flow_through_potential = Q.*params.Q0 .*(psi + dNdx);
    subplot(3,2,6);
    plot(x_grid_hydro./1000,melt_term,'DisplayName', strcat(I_F,': ','Melt'));
    hold on;
    plot(x_grid_hydro./1000,flow_through_potential,'DisplayName', strcat(I_F,': ','Flow through potential'));
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$mL = Q(\psi + \frac{\partial N}{\partial x})$','Interpreter','latex')
end 


%% Helper functions
function b = bed(x,params)
    xsill = x>params.sill_min & x<params.sill_max;
    xdsill = x>=params.sill_max;
    sill_length = params.sill_max-params.sill_min;
    
    b = params.b0 + params.bx.*x;
    
    b(xsill) = params.b0 + (params.bx*params.sill_min) + params.sill_slope.*(x(xsill)-params.sill_min);

    b(xdsill) = params.b0 + (params.bx*params.sill_min) + params.sill_slope.*sill_length + ...
            params.bx*(x(xdsill)-params.sill_max);

end