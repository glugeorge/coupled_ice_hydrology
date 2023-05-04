figure();
t = tiledlayout(3,3,'TileSpacing','Compact');
load slab_oneway_1.mat;
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 
xg = params.xg;
h = params.h;
u = params.u;
Q = Qs(end,:)';
N = Ns(end,:)';
S = Ss(end,:)';
x_grid_hydro = sigma_h*xg*params.x0;
h_interp = interp1(sigma_elem,h.*params.h0,sigma_h,"linear","extrap");
u_interp = interp1(sigma,u.*params.u0,sigma_h,"linear","extrap");
b = -bed(x_grid_hydro,params);
psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);

nexttile;
plot(x_grid_hydro./1000,params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro),'-g','LineWidth',2,'DisplayName','$\rho_wg\frac{\partial b}{\partial x}$');
hold on; 
plot(x_grid_hydro./1000,gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro),'-b','LineWidth',2,'DisplayName','$\rho_ig\frac{\partial h}{\partial x}$')
plot(x_grid_hydro./1000,psi,'k','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');ylim([-1 11]);
legend('Interpreter','latex');
title('$\psi = \rho_wg\frac{\partial b}{\partial x} - \rho_ig\frac{\partial h}{\partial x}$','Interpreter','latex');
ylabel("Slab",'Interpreter','latex');

nexttile;
m_over_rho_i = params.m0.*abs(Q).^3./(params.rho_i.*S.^(8/3)); % from m'=|Q'|^3/S'^8/3
KSN3 = params.K0.*S.*params.S0.*(N.*params.N0).^3;
u_advect = u_interp.*gradient(S.*params.S0)./gradient(x_grid_hydro);
plot(x_grid_hydro./1000,m_over_rho_i,'r','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{m}{\rho_i}$');
hold on;
plot(x_grid_hydro./1000,KSN3,'--k','LineWidth',1,'MarkerSize',8,'DisplayName', '$K_0SN^3$');
plot(x_grid_hydro./1000,u_advect,':k','LineWidth',2,'MarkerSize',8,'DisplayName', '$u\frac{\partial S}{\partial x}$');
legend('Interpreter','latex');
title('$\frac{\partial S}{\partial t_h} = \frac{m}{\rho_i}-K_0SN^3 - u\frac{\partial S}{\partial x}$',Interpreter='latex')

nexttile;
dNdx = gradient(N.*params.N0)./gradient(x_grid_hydro);
f_Q_S = params.f*params.rho_w*params.g.*Q.*params.Q0.*abs(Q.*params.Q0)./(S.*params.S0).^(8/3);
h_interp = interp1(sigma_elem,h,sigma_h,'linear','extrap');
b = -bed(x_grid_hydro,params);
psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);
plot(x_grid_hydro./1000,psi,'k','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');
hold on;
plot(x_grid_hydro./1000,dNdx,'-.','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{\partial N}{\partial x}$');
plot(x_grid_hydro./1000,f_Q_S,':','LineWidth',2,'MarkerSize',8,'DisplayName', '$f\rho_wg\frac{Q|Q|}{S^{8/3}}$');ylim([-10 11])
legend('Interpreter','latex');
title('$\psi+ \frac{\partial N}{\partial x} = f\rho_wg\frac{Q|Q|}{S^{8/3}}$','Interpreter','latex')

% Quadratic
load quad_oneway_1.mat;
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 
xg = params.xg;
h = params.h;
u = params.u;
Q = Qs(end,:)';
N = Ns(end,:)';
S = Ss(end,:)';
x_grid_hydro = sigma_h*xg*params.x0;
h_interp = interp1(sigma_elem,h.*params.h0,sigma_h,"linear","extrap");
u_interp = interp1(sigma,u.*params.u0,sigma_h,"linear","extrap");
b = -bed(x_grid_hydro,params);
psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);

nexttile;
plot(x_grid_hydro./1000,params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro),'-g','LineWidth',2,'MarkerSize',8,'DisplayName','$\rho_wg\frac{\partial b}{\partial x}$');
hold on; 
plot(x_grid_hydro./1000,gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro),'-b','LineWidth',2,'MarkerSize',8,'DisplayName','$\rho_ig\frac{\partial h}{\partial x}$')
plot(x_grid_hydro./1000,psi,'k','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');
ylabel('Quadratic','Interpreter','latex');
nexttile;
m_over_rho_i = params.m0.*abs(Q).^3./(params.rho_i.*S.^(8/3)); % from m'=|Q'|^3/S'^8/3
KSN3 = params.K0.*S.*params.S0.*(N.*params.N0).^3;
u_advect = u_interp.*gradient(S.*params.S0)./gradient(x_grid_hydro);
plot(x_grid_hydro./1000,m_over_rho_i,'r','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{m}{\rho_i}$');
hold on;
plot(x_grid_hydro./1000,KSN3,'--k','LineWidth',1,'MarkerSize',8,'DisplayName', '$K_0SN^3$');
plot(x_grid_hydro./1000,u_advect,':k','LineWidth',2,'MarkerSize',8,'DisplayName', '$u\frac{\partial S}{\partial x}$');

nexttile;
dNdx = gradient(N.*params.N0)./gradient(x_grid_hydro);
f_Q_S = params.f*params.rho_w*params.g.*Q.*params.Q0.*abs(Q.*params.Q0)./(S.*params.S0).^(8/3);
h_interp = interp1(sigma_elem,h,sigma_h,'linear','extrap');
b = -bed(x_grid_hydro,params);
psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);
plot(x_grid_hydro./1000,psi,'k','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');
hold on;
plot(x_grid_hydro./1000,dNdx,'-.','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{\partial N}{\partial x}$');
plot(x_grid_hydro./1000,f_Q_S,':','LineWidth',2,'MarkerSize',8,'DisplayName', '$f\rho_wg\frac{Q|Q|}{S^{8/3}}$');

% Coupled
load steady_state_combined.mat
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 
x_grid_hydro = sigma_h*xg*params.x0;
h_interp = interp1(sigma_elem,h.*params.h0,sigma_h,"linear","extrap");
u_interp = interp1(sigma,u.*params.u0,sigma_h,"linear","extrap");
b = -bed(x_grid_hydro,params);
psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);

nexttile;
plot(x_grid_hydro./1000,params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro),'-g','LineWidth',2,'MarkerSize',8,'DisplayName','$\rho_wg\frac{\partial b}{\partial x}$');
hold on; 
plot(x_grid_hydro./1000,gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro),'-b','LineWidth',2,'MarkerSize',8,'DisplayName','$\rho_ig\frac{\partial h}{\partial x}$')
plot(x_grid_hydro./1000,psi,'k','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');
ylabel('Coupled','Interpreter','latex')

nexttile;
m_over_rho_i = params.m0.*abs(Q).^3./(params.rho_i.*S.^(8/3)); % from m'=|Q'|^3/S'^8/3
KSN3 = params.K0.*S.*params.S0.*(N.*params.N0).^3;
u_advect = u_interp.*gradient(S.*params.S0)./gradient(x_grid_hydro);
plot(x_grid_hydro./1000,m_over_rho_i,'r','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{m}{\rho_i}$');
hold on;
plot(x_grid_hydro./1000,KSN3,'--k','LineWidth',1,'MarkerSize',8,'DisplayName', '$K_0SN^3$');
plot(x_grid_hydro./1000,u_advect,':k','LineWidth',2,'MarkerSize',8,'DisplayName', '$u\frac{\partial S}{\partial x}$');

nexttile;
dNdx = gradient(N.*params.N0)./gradient(x_grid_hydro);
f_Q_S = params.f*params.rho_w*params.g.*Q.*params.Q0.*abs(Q.*params.Q0)./(S.*params.S0).^(8/3);
h_interp = interp1(sigma_elem,h,sigma_h,'linear','extrap');
b = -bed(x_grid_hydro,params);
psi = params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro);
plot(x_grid_hydro./1000,psi,'k','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');
hold on;
plot(x_grid_hydro./1000,dNdx,'-.','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{\partial N}{\partial x}$');
plot(x_grid_hydro./1000,f_Q_S,':','LineWidth',2,'MarkerSize',8,'DisplayName', '$f\rho_wg\frac{Q|Q|}{S^{8/3}}$');
xlabel(t,'distance from divide, \emph{x} (km)','Interpreter','latex');

figure()
x_grid_ice = sigma*xg*params.x0;
h_interp = interp1(sigma_elem,h.*params.h0,sigma,"linear","extrap");
N_interp = interp1(sigma_h,N,sigma,"linear","extrap");
shear_stress = params.C.*N_interp(2:end)*params.N0.*(u(2:end).*params.u0./(u(2:end).*params.u0+params.As*(params.C*N_interp(2:end)*params.N0).^params.n)).^(1/params.n); 
b = -bed(x_grid_ice,params);
driving_stress = params.rho_i*params.g*h_interp(2:end).*gradient(h_interp(2:end)-b(2:end))./gradient(x_grid_ice(2:end));
shear_and_driving = shear_stress + driving_stress;

dudx = gradient(u(2:end).*params.u0)./gradient(x_grid_ice(2:end));
long_stress = gradient(2*params.A^(-1/params.n)*h_interp(2:end).*abs(dudx).^(1/params.n -1).*dudx)./gradient(x_grid_ice(2:end));
plot(x_grid_ice(2:end)./1000,shear_stress,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Basal Shear');
hold on;
plot(x_grid_ice(2:end)./1000,driving_stress,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Driving');
plot(x_grid_ice(2:end)./1000,long_stress,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Longitudinal');
legend('Interpreter','latex');
title('Stress terms from Equation 11','Interpreter','latex')

xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
ylabel('stress (Pa)','Interpreter','latex')
