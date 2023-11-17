fig = figure(Name='Figure 5',NumberTitle='off',Position=[159 93 1200,675]);
t = tiledlayout(3,3,'TileSpacing','Compact');
set(groot,'defaultAxesTickLabelInterpreter','latex');  

% Budd coupled
load budd_steady_state.mat;
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 

x_grid_hydro = sigma_h*xg*params.x0;
h_interp = interp1(sigma_elem,h,sigma_h,"linear","extrap");
u_interp = interp1(sigma,u,sigma_h,"linear","extrap");
b = -bed(x_grid_hydro,params);

nexttile;
melt = abs(Q).^3./(S.^(8/3)); 
creep = S.*N.^3;
advect = params.beta.*u_interp.*gradient(S)./gradient(sigma_h.*xg);
plot(sigma_h,melt,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{|Q|^3}{S^{8/3}}$');
hold on;
plot(sigma_h,creep,'--','LineWidth',1,'MarkerSize',8,'DisplayName', '$SN^3$');
plot(sigma_h,advect,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\beta u\frac{\partial S}{\partial x}$');
legend('Interpreter','latex','Location','northwest');
title('$\frac{\partial S}{\partial t_h} = \frac{|Q|^3}{S^{8/3}}-SN^3 - \beta u\frac{\partial S}{\partial x}$',Interpreter='latex')
ylabel('S1.B','Interpreter','latex');
%xlim([0.95 1]);ylim([0 0.6]);

nexttile;
delta_dNdx = params.delta.*gradient(N)./gradient(sigma_h.*xg);
Q_S = Q.*abs(Q)./(S).^(8/3);
psi = (params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro))./params.psi0;
plot(sigma_h,psi,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');
hold on;
plot(sigma_h,delta_dNdx,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\delta\frac{\partial N}{\partial x}$');
plot(sigma_h,Q_S,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{Q|Q|}{S^{8/3}}$');
legend('Interpreter','latex','Location','northwest');
title('$\psi+ \delta\frac{\partial N}{\partial x} = \frac{Q|Q|}{S^{8/3}}$','Interpreter','latex')

nexttile;
tau = N.*u_interp.^(1/params.n); 
b = -bed(x_grid_hydro,params)./params.h0;
driving_stress = h_interp.*gradient(h_interp-b)./gradient(sigma_h.*xg);

dudx = gradient(u_interp)./gradient(sigma_h.*xg);
long_stress = params.alpha.*gradient(h_interp.*abs(dudx).^(1/params.n -1).*dudx)./gradient(sigma_h.*xg);
plot(sigma_h,tau,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Basal Shear');
hold on;
plot(sigma_h,driving_stress,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Driving');
plot(sigma_h,long_stress,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Longitudinal');
legend('Interpreter','latex','Location','northwest');
title('Stress terms in Eq. (16)','Interpreter','latex');
% Coulomb coupled
load coulomb_steady_state.mat;
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 

x_grid_hydro = sigma_h*xg*params.x0;
h_interp = interp1(sigma_elem,h,sigma_h,"linear","extrap");
u_interp = interp1(sigma,u,sigma_h,"linear","extrap");
b = -bed(x_grid_hydro,params);
legend('Interpreter','latex','Location','northwest');


nexttile;
melt = abs(Q).^3./(S.^(8/3)); 
creep = S.*N.^3;
advect = params.beta.*u_interp.*gradient(S)./gradient(sigma_h.*xg);
plot(sigma_h,melt,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{|Q|^3}{S^{8/3}}$');
hold on;
plot(sigma_h,creep,'--','LineWidth',1,'MarkerSize',8,'DisplayName', '$SN^3$');
plot(sigma_h,advect,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\beta u\frac{\partial S}{\partial x}$');
ylabel('S1.C','Interpreter','latex');
legend('Interpreter','latex','Location','northwest');

%xlim([0.95 1]);ylim([0 0.6]);

nexttile;
delta_dNdx = params.delta.*gradient(N)./gradient(sigma_h.*xg);
Q_S = Q.*abs(Q)./(S).^(8/3);
psi = (params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro))./params.psi0;
plot(sigma_h,psi,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');
hold on;
plot(sigma_h,delta_dNdx,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{\partial N}{\partial x}$');
plot(sigma_h,Q_S,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$f\rho_wg\frac{Q|Q|}{S^{8/3}}$');
legend('Interpreter','latex','Location','northwest');

nexttile;
tau = params.gamma.*N.*(u_interp./(u_interp+N.^params.n)).^(1/params.n); 
b = -bed(x_grid_hydro,params)./params.h0;
driving_stress = h_interp.*gradient(h_interp-b)./gradient(sigma_h.*xg);

dudx = gradient(u_interp)./gradient(sigma_h.*xg);
long_stress = params.alpha.*gradient(h_interp.*abs(dudx).^(1/params.n -1).*dudx)./gradient(sigma_h.*xg);
plot(sigma_h,tau,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Basal Shear');
hold on;
plot(sigma_h,driving_stress,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Driving');
plot(sigma_h,long_stress,'LineWidth',1,'MarkerSize',8,'DisplayName', 'Longitudinal');
legend('Interpreter','latex','Location','northwest');

% Quadratic
load quad_oneway_1.mat;
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 

x_grid_hydro = sigma_h*params.xg*params.x0;
h_interp = interp1(sigma_elem,params.h,sigma_h,"linear","extrap");
u_interp = interp1(sigma,params.u,sigma_h,"linear","extrap");
Q = Qs';
S = Ss';
N = Ns';
b = -bed(x_grid_hydro,params);
legend('Interpreter','latex','Location','northwest');

nexttile;
melt = abs(Q).^3./(S.^(8/3)); 
creep = S.*N.^3;
advect = params.beta.*u_interp.*gradient(S)./gradient(sigma_h.*xg);
plot(sigma_h,melt,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{|Q|^3}{S^{8/3}}$');
hold on;
plot(sigma_h,creep,'--','LineWidth',1,'MarkerSize',8,'DisplayName', '$SN^3$');
plot(sigma_h,advect,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\beta u\frac{\partial S}{\partial x}$');
ylabel('S2','Interpreter','latex');
%xlim([0.95 1]);ylim([0 0.6]);
legend('Interpreter','latex','Location','northwest');

nexttile;
delta_dNdx = params.delta.*gradient(N)./gradient(sigma_h.*xg);
Q_S = Q.*abs(Q)./(S).^(8/3);
psi = (params.rho_w*params.g*gradient(b)./gradient(x_grid_hydro) - gradient(params.rho_i*params.g*h_interp*params.h0)./gradient(x_grid_hydro))./params.psi0;
plot(sigma_h,psi,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\psi$');
hold on;
plot(sigma_h,delta_dNdx,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$\frac{\partial N}{\partial x}$');
plot(sigma_h,Q_S,'-','LineWidth',1,'MarkerSize',8,'DisplayName', '$f\rho_wg\frac{Q|Q|}{S^{8/3}}$');
legend('Interpreter','latex','Location','northwest');

nexttile;
xlabel(t,'Normalized distance from divide, $\sigma$','Interpreter','latex')
fontsize(fig, 14, "points")
