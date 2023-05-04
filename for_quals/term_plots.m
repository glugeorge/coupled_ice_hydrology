close all; clear;
load schoof_retreat_highres.mat
QNShuxg = results.init_cond;
params = results.params;
plot_terms(QNShuxg,params,'Initial')
%  for t=1:10:100
%      plot_vars_new(results,t,params)
%  end

%QNShuxg = results.steady_state;
% params = results.params;
% plot_terms(QNShuxg,params,'Final')
%figure();
%load C_0.2_A_4.9_c_nosplit_highres.mat;
%QNShuxg = results.init_cond;
%params = results.params;
%plot_vars_new(results,1,params)

%plot_terms(QNShuxg,params,'Initial')

function plot_terms(QNShuxg,params,I_F)
    sigma_elem = params.sigma_elem;
    sigma_h = params.sigma_h;
    sigma = params.sigma; 
    N1 = params.N1;
%     Q = QNShuxg(1:params.Nx);
%     N = QNShuxg(params.Nx+1:2*params.Nx);
%     S = QNShuxg(2*params.Nx+1:3*params.Nx);
%     h = QNShuxg(3*params.Nx+1:4*params.Nx);
%     u = QNShuxg(4*params.Nx+1:5*params.Nx);
%     xg = QNShuxg(5*params.Nx+1);
    Q = QNShuxg(1:params.Nh);
    N = QNShuxg(params.Nh+1:2*params.Nh);
    S = QNShuxg(2*params.Nh+1:3*params.Nh);
    h = QNShuxg(params.ice_start+1:params.ice_start+ params.Nx);
    u = QNShuxg(params.ice_start + params.Nx+1:params.ice_start+2*params.Nx);
    xg = QNShuxg(params.ice_start+2*params.Nx+1);
    
    % do ice plots on sigma_elem grid
    x_grid_ice = sigma*xg*params.x0;
    diff_ice = diff(x_grid_ice);
    d_coarse = diff_ice(1);
    d_fine = diff_ice(N1);
    figure();
    %% Plot 1 - ice mass conservation
    % Remember that since these are all in steady state, d/dt goes to zero
    accum = params.accum; % normalization done in flowline equations
    h_interp = interp1(sigma_elem,h.*params.h0,sigma,"linear","extrap");
    dhu_dx = zeros(size(x_grid_ice));
    dhu_dx(1:N1) = gradient((h_interp(1:N1).*u(1:N1)).*params.u0)./d_coarse;
    dhu_dx(N1+1:end) = gradient((h_interp(N1+1:end).*u(N1+1:end)).*params.u0)./d_fine;
    subplot(3,2,1);
    plot(x_grid_ice./1000,accum.*ones(size(x_grid_ice)),'DisplayName', strcat(I_F,': ','Accumulation'));
    hold on;
    plot(x_grid_ice./1000,dhu_dx,'--','DisplayName', strcat(I_F,': ','$\frac{\partial (hu)}{\partial x}$'));
    legend('Interpreter','latex');
    xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    title('$\frac{\partial h_g}{\partial t} + \frac{\partial (h_gu)}{\partial x}= a$','Interpreter','latex')
    %% Plot 2: ice momentum conservation
    N_interp = interp1(sigma_h,N,sigma,"linear","extrap");
    shear_stress = params.C.*N_interp(2:end)*params.N0.*(u(2:end).*params.u0./(u(2:end).*params.u0+params.As*(params.C*N_interp(2:end)*params.N0).^params.n)).^(1/params.n); 
    b = -bed_schoof(x_grid_ice,params);
    driving_stress = params.rho_i*params.g*h_interp(2:end).*gradient(h_interp(2:end)-b(2:end))./gradient(x_grid_ice(2:end));
    shear_and_driving = shear_stress + driving_stress;
    
    dudx = gradient(u(2:end).*params.u0)./gradient(x_grid_ice(2:end));
    long_stress = gradient(2*params.A^(-1/params.n)*h_interp(2:end).*abs(dudx).^(1/params.n -1).*dudx)./gradient(x_grid_ice(2:end));
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

function plot_vars(QNShuxg,params)
    sigma_elem = params.sigma_elem;
    sigma_h = params.sigma_h;
    sigma = params.sigma; 
    Q = QNShuxg(1:params.Nx);
    N = QNShuxg(params.Nx+1:2*params.Nx);
    S = QNShuxg(2*params.Nx+1:3*params.Nx);
    h = QNShuxg(3*params.Nx+1:4*params.Nx);
    u = QNShuxg(4*params.Nx+1:5*params.Nx);
    xg = QNShuxg(5*params.Nx+1);
    
    figure(); 
    ax1 = subplot(5,1,1); 
    plot(sigma_h.*xg.*params.x0./1000,Q); ylabel('Q');title('Nondimensionalized steady state variables');
    ax2 =subplot(5,1,2);
    plot(sigma_h.*xg.*params.x0./1000,N);ylabel('N');
    ax3 = subplot(5,1,3);
    plot(sigma_h.*xg.*params.x0./1000,S);ylabel('S');
    ax4 = subplot(5,1,4);
    plot(sigma_elem.*xg.*params.x0./1000,h);ylabel('h');
    ax5 = subplot(5,1,5);
    plot(sigma.*xg.*params.x0./1000,u);ylabel('u');xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    linkaxes([ax1,ax2,ax3,ax4,ax5],'x')

    

    
end

function plot_vars_new(results,time,params)
    sigma_elem = params.sigma_elem;
    sigma_h = params.sigma_h;
    sigma = params.sigma; 
    Q = results.Qs(:,time);
    N = results.Ns(:,time);
    S = results.Ss(:,time);
    h = results.hs(:,time);
    u = results.us(:,time);
    xg = results.xgs(time);
        
    ax1 = subplot(5,1,1); 
    plot(sigma_h.*xg.*params.x0./1000,Q); ylabel('Q');title('Nondimensionalized steady state variables');
    hold on;
    if mod(time-1,20) == 0
    text(max(sigma_h.*xg.*params.x0./1000),max(Q),['t=',num2str(time-1)])
    end
    ax2 =subplot(5,1,2);
    plot(sigma_h.*xg.*params.x0./1000,N);ylabel('N');hold on;
    ax3 = subplot(5,1,3);
    plot(sigma_h.*xg.*params.x0./1000,S);ylabel('S');hold on;
    ax4 = subplot(5,1,4);
    plot(sigma_elem.*xg.*params.x0./1000,h);ylabel('h');hold on;
    if time==81
        plot(linspace(0,1)*1500,bed_schoof(linspace(0,1).*1500e3,params)./params.h0,'-k'); 
    end
    ax5 = subplot(5,1,5);
    plot(sigma.*xg.*params.x0./1000,u);hold on;ylabel('u');xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
    xlim([0,1500]);

    

    
end