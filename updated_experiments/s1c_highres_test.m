clear; close all;
test1 = 1; %vary between 0.5 and 1.5
test2 = 1; % vary between 0.5 and 2.5
test3 = 0.3; % vary between 0.3 and 0.6

%% Define parameter ranges
A_vals = logspace(log10(3.9e-26),log10(4.2e-25),4);
M_vals = logspace(-5,-3,4);
a_vals = linspace(0.1,0.5,4);

% C values
C_vals = linspace(0.1,0.5,4);

% Base vector
A = median(A_vals);
M = median(M_vals);
accum = median(a_vals);
C = median(C_vals);

%% Bed parameters
params.b0 = -100;           %bed topo at x=0
params.bx = -1e-3;          %linear bed slope

params.sill_min = 2000e3;   %sill min x position
params.sill_max = 2100e3;   %sill max x position
params.sill_slope = 1e-3;   %slope of sill
  
%% Physical parameters
params.A = A; 
params.n = 3;
params.rho_i = 917;
params.rho_w = 1028;
params.g = 9.81;
params.C = C; % base 0.2
params.As = 2.26e-21; % Calculated 
params.f = 0.07; % From Kingslake thesis
params.K0 = 10^-24; % From Kingslake thesis  
params.L = 3.3e5; % Kingslake thesis
params.year = 3600*24*365;
%% Scaling params (coupled model equations solved in non-dim form)
params.x0 = 100*10^3;
params.h0 = 1000;
params.Q0 = 1;

params.psi0 = params.rho_w*params.g*params.h0/params.x0;
params.M0 = params.Q0/params.x0;
params.m0 = params.Q0*params.psi0/params.L;
params.eps_r = params.m0*params.x0/(params.rho_i*params.Q0);
params.S0 = (params.f*params.rho_w*params.g*params.Q0^2/params.psi0)^(3/8);
params.th0 = params.rho_i*params.S0/params.m0;
params.N0 = (params.K0*params.th0)^(-1/3);
params.delta = params.N0/(params.x0*params.psi0);
%params.u0 = (rho_i*g*params.h0^2/(C*params.N0*params.x0))^n;
params.u0 = params.As*(params.C*params.N0)^params.n;
params.t0 = params.x0/params.u0;
params.a0 = params.h0/params.t0;
params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));
%gamma = As*(C*N0)^n;
params.gamma = (params.C*params.N0*params.x0)/(params.rho_i*params.g*params.h0^2);
params.beta = params.th0/params.t0;
params.r = params.rho_i/params.rho_w;

params.transient = 0;
%% Establish timings
params.year = 3600*24*365;  %number of seconds in a year
params.Nt =50;                    %number of time steps - normally 150
params.end_year = 1000; %normally 7500

params.dt = params.end_year*params.year/params.Nt;

%% Determine at what points there is coupling
% 1 - coupling on, 0 - coupling off
params.hydro_u_from_ice_u = 1;
params.hydro_psi_from_ice_h = 1;
params.ice_N_from_hydro = 1;
params.Cf = 1;

% Grid parameters - ice sheet
params.Nx = 3000;                    %number of grid points - 200
params.Nh = 3000;
fine_grid = linspace(0,1,2*params.Nx);
params.sigma = fine_grid(2:2:end)';    %grid points on velocity (includes GL, not ice divide)
params.sigma_elem = fine_grid(1:2:end)';
params.dsigma = diff(params.sigma_elem); %grid spacing

%% Grid parameters - hydro
params.sigma_h = linspace(0,1,params.Nh)';
params.dsigma_h = diff(params.sigma_h); %grid spacing


%% Initial "steady state" conditions
params.shear_scale = 1;
Q = ones(params.Nh,1);
N = test1*ones(params.Nh,1);
S = ones(params.Nh,1); 
%Q = interp1(sigma_h_old,Q_old,params.sigma_h,"linear","extrap");
%N = interp1(sigma_h_old,N_old,params.sigma_h,"linear","extrap");
%S = interp1(sigma_h_old,S_old,params.sigma_h,"linear","extrap");

params.S_old = S;
params.M = M./params.M0; % zero when using schoof bed
params.N_terminus = 0;
params.accum = accum./params.year;
xg = 200e3/params.x0;
%xg = xg_old;
%hf = (-bed_flat(xg.*params.x0,params)/params.h0)/params.r;
h =  test2*ones(length(params.sigma_elem),1);%1 - (1-hf).*params.sigma;
%h = interp1(sigma_elem_old,h_old,params.sigma_elem,"linear","extrap");
u = test3*(params.sigma_elem.^(1/3)) + 1e-3; % 0.1 for C = 0.5, 0.3 for C = 0.1-0.4
%u = interp1(sigma_old,u_old,params.sigma,"linear","extrap");

params.Q_in = 0.001/params.Q0;

params.h_old = h;
params.xg_old = xg;
params.ice_start = 3*params.Nh;

sig_old = params.sigma;
sige_old = params.sigma_elem;
QNShuxg0 = [Q;N;S;h;u;xg];

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
flf = @(QNShuxg) combined_hydro_ice_eqns_highres(QNShuxg,params);

[QNShuxg_init,F,exitflag,output,JAC] = fsolve(flf,QNShuxg0,options);

Q = QNShuxg_init(1:params.Nh);
N = QNShuxg_init(params.Nh+1:2*params.Nh);
S = QNShuxg_init(2*params.Nh+1:3*params.Nh);
h = QNShuxg_init(params.ice_start+1:params.ice_start+ params.Nx);
u = QNShuxg_init(params.ice_start + params.Nx+1:params.ice_start+2*params.Nx);
xg = QNShuxg_init(params.ice_start+2*params.Nx+1);
hf = (-bed_flat(xg.*params.x0,params)/params.h0)/(params.r);
%% Save
save("S1C_highres.mat");

%% Plot terms
%plot_terms(QNShuxg_init,params);
figure()
ax1 = subplot(5,1,1); 
plot(params.sigma_h.*xg.*params.x0./1000,Q); ylabel('Q');title('Nondimensionalized steady state variables');

ax2 =subplot(5,1,2);
plot(params.sigma_h.*xg.*params.x0./1000,N);ylabel('N');hold on;
ax3 = subplot(5,1,3);
plot(params.sigma_h.*xg.*params.x0./1000,S);ylabel('S');hold on;
ax4 = subplot(5,1,4);
b_h = -bed_flat(xg.*params.sigma_elem.*params.x0,params)./params.h0;
    
plot(params.sigma_elem.*xg.*params.x0./1000,h-b_h);ylabel('h');hold on;plot(params.sigma_elem.*xg.*params.x0./1000,-b_h,'k')
ax5 = subplot(5,1,5);
plot(params.sigma.*xg.*params.x0./1000,u);hold on;ylabel('u');xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
linkaxes([ax1,ax2,ax3,ax4,ax5],'x')



function F = combined_hydro_ice_eqns_highres(QNShuxg,params)
    % only for steady state solutions
    % unpack variables
    M = params.M;
    Q_in = params.Q_in;
    Nx = params.Nx;
    Nh = params.Nh;
    ice_start = params.ice_start;
    Q = QNShuxg(1:Nh);
    N = QNShuxg(Nh+1:2*Nh);
    S = QNShuxg(2*Nh+1:3*Nh);
    h = QNShuxg(ice_start+1:ice_start + Nx);
    u = QNShuxg(ice_start + Nx+1:ice_start + 2*Nx);
    xg = QNShuxg(ice_start + 2*Nx+1);
    hf = (-bed_flat(xg.*params.x0,params)/params.h0)/(params.r);

    %grid params unpack
    dt = params.dt/params.t0;
    dth= params.dt/params.th0;
    ds = params.dsigma;
    Nx = params.Nx;
    sigma = params.sigma;
    sigma_elem = params.sigma_elem;
    b = -bed_flat(xg.*sigma.*params.x0,params)/params.h0;
    dbdx = gradient(bed_flat(xg.*params.sigma_h.*params.x0,params))./gradient(params.sigma_h.*params.x0*xg);
    Fh = zeros(Nx,1);
    Fu = zeros(Nx,1);
    
    %physical params unpack
    m     = 1/params.n;
    nglen = params.n;
    accum = params.accum;
    a = accum/params.a0;
    gamma = params.gamma;
    alpha = params.alpha;
    ss = params.transient;
       
    %previous time step unpack
    h_old = params.h_old;
    xg_old = params.xg_old;
    
    % Assign values based off coupling parameters
    if params.hydro_u_from_ice_u
        u_ice_interp = interp1(sigma,u,params.sigma_h,"linear","extrap");
    else
        u_ice_interp = 1000.*ones(params.Nx,1)./(params.year*params.u0);
    end

    if params.hydro_psi_from_ice_h
        %gradient_h = zeros(size(sigma_elem));
        %gradient_h(1:N1) = gradient(h(1:N1).*params.h0)./gradient(sigma_elem(1:N1).*params.x0*xg);
        %gradient_h(N1+1:end) = gradient(h(N1+1:end).*params.h0)./gradient(sigma_elem(N1+1:end).*params.x0*xg);
        %gradient_h_interp = interp1(sigma_elem,gradient_h,params.sigma_h,'linear','extrap');
        %params.psi = (-params.rho_w*params.g.*dbdx-params.rho_i*params.g.*gradient_h_interp)./params.psi0;
        h_interp = interp1(sigma_elem,h.*params.h0,params.sigma_h,'linear','extrap');
        params.psi = (-params.rho_w*params.g.*dbdx-params.rho_i*params.g.*gradient(h_interp)./gradient(params.sigma_h.*params.x0*xg))./params.psi0;
    else
        params.psi = 1*(1-3*exp(-20.*params.sigma_h));
    end

    if params.ice_N_from_hydro
        N_ice = interp1(params.sigma_h,N,sigma);
    else
        N_ice = 100000*ones(size(h_old))./params.N0;
    end

    % Q
    fq(1) = Q(1) - Q_in; % Boundary condition

    fq(2:Nh-1) = (params.eps_r*(params.r-1))*abs(Q(2:Nh-1)).^3./S(2:Nh-1).^(8/3) +...
                        params.eps_r.*S(2:Nh-1).*N(2:Nh-1).^3 + M + ...
                        params.eps_r.*params.beta.*u_ice_interp(2:Nh-1).*(S(3:Nh)-S(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1)) - ...
                        (Q(3:Nh)-Q(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1));
    fq(Nh) = (params.eps_r*(params.r-1))*abs(Q(Nh)).^3./S(Nh).^(8/3) +...
                        params.eps_r*S(Nh).*N(Nh).^3 + M + ...
                        params.eps_r.*params.beta.*u_ice_interp(Nh).*(S(Nh)-S(Nh-1))./(xg*params.dsigma_h(Nh-1)) - ...
                        (Q(Nh)-Q(Nh-1))./(xg*params.dsigma_h(Nh-1)); % one sided difference instead
    % N 
    fn(1) = Q(1).*abs(Q(1))./(S(1).^(8/3)) - params.psi(1) - ...
                        params.delta*(N(2)-N(1))./(xg*params.dsigma_h(1)); % use 1 sided difference instead of symmetry argument
    fn(2:Nh-1) =  Q(2:Nh-1).*abs(Q(2:Nh-1))./(S(2:Nh-1).^(8/3)) - params.psi(2:Nh-1) - ...
                        params.delta*(N(3:Nh)-N(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1));
    fn(Nh) = N(Nh) - params.N_terminus; % Boundary condition
    % S
    fs(1) = abs(Q(1)).^3./(S(1).^(8/3)) - ... 
                        S(1).*N(1).^3 + ...
                        0; % one sided difference
    fs(2:Nh-1)= abs(Q(2:Nh-1)).^3./(S(2:Nh-1).^(8/3)) - ... 
                        S(2:Nh-1).*N(2:Nh-1).^3 + ...
                        (0 - params.beta.*u_ice_interp(2:Nh-1)).*(S(3:Nh)-S(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1)) - ...
                       0; 
    fs(Nh)= abs(Q(Nh)).^3./(S(Nh).^(8/3)) - ... 
                        S(Nh).*N(Nh).^3 + ...
                        (0 - params.beta.*u_ice_interp(Nh)).*(S(Nh)-S(Nh-1))./(xg*params.dsigma_h(Nh-1)) - ...
                        0; % one sided difference

    % ice sheet equations
    Fh(1)      = 0 + (2.*h(1).*u(1))./(ds(1).*xg) - a;
    Fh(2)      = 0 -...
                   (h(2).*(u(2)+u(1)))./(2*xg.*ds(2)) -...
                            a;
    Fh(3:Nx-1) = 0 -...
                   (h(3:Nx-1).*(u(3:Nx-1)+u(2:Nx-2)) - h(2:Nx-2).*(u(2:Nx-2)+u(1:Nx-3)))./(2*xg.*ds(3:Nx-1)) -...
                            a;
                            
    Fh(Nx)     = 0 -...
                        (h(Nx).*(u(Nx)+u(Nx-1)) - h(Nx-1).*(u(Nx-1)+u(Nx-2)))./(2*xg.*ds(Nx-1)) -...
                            a;

	%velocity
    Fu(1)      = (alpha).*(1./(xg.*ds(1)).^((1/nglen)+1)).*...
                 (h(2).*(u(2)-u(1)).*abs(u(2)-u(1)).^((1/nglen)-1) -...
                  h(1).*(2*u(1)).*abs(2*u(1)).^((1/nglen)-1)) -...
                  gamma.*params.shear_scale.*N_ice(1).* (u(1)./(u(1) + N_ice(1).^nglen)).^m -...
                  0.5.*(h(1)+h(2)).*(h(2)-b(2)-h(1)+b(1))./(xg.*ds(1));
    Fu(2:Nx-1) = (alpha).*(1./(xg.*ds(2:Nx-1)).^((1/nglen)+1)).*...
                 (h(3:Nx).*(u(3:Nx)-u(2:Nx-1)).*abs(u(3:Nx)-u(2:Nx-1)).^((1/nglen)-1) -...
                  h(2:Nx-1).*(u(2:Nx-1)-u(1:Nx-2)).*abs(u(2:Nx-1)-u(1:Nx-2)).^((1/nglen)-1)) -...
                  gamma.*params.shear_scale.*N_ice(2:Nx-1).*(u(2:Nx-1)./(u(2:Nx-1) + N_ice(2:Nx-1).^nglen)).^m -...
                  0.5.*(h(2:Nx-1)+h(3:Nx)).*(h(3:Nx)-b(3:Nx)-h(2:Nx-1)+b(2:Nx-1))./(xg.*ds(2:Nx-1));
    Fu(Nx)     = alpha.*(1./(xg.*ds(Nx-1)).^(1/nglen)).*...
                 (abs(u(Nx)-u(Nx-1)).^((1/nglen)-1)).*(u(Nx)-u(Nx-1)) - 0.5*(1-params.r)*hf;
             
    Fxg        = 3*h(Nx) - h(Nx-1) - 2*hf; 

    F = [fq';fn';fs';Fh;Fu;Fxg];
end
%% helper functions
function plot_terms(QNShuxg,params)
    sigma_elem = params.sigma_elem;
    sigma_h = params.sigma_h;
    sigma = params.sigma; 
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
