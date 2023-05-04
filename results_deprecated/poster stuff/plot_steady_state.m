% First plot runs without coupling 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

figure("Color",'white','Position',[440,136,887,661])
title('Steady state solutions for ice flow')
ax1 = subplot('Position',[.1 .81 .8 .17]);
% First plot ice stuff
% Plot Ice thickness changes with coupling
load SSA_simple_result.mat
plot(params.sigma_elem*xg*params.xscale/1000,h*params.hscale,'DisplayName','Power Law','LineWidth',3,'Color','r');
hold on;
load SSA_const_N_result.mat
plot(params.sigma_elem*xg/1000,h*params.h0,'DisplayName','Generalized Coulomb Law, uniform N','LineWidth',2,'Color','k','LineStyle','--');
load SSA_sea_level_result.mat
plot(params.sigma_elem*xg/1000,h*params.h0,'DisplayName','Generalized Coulomb Law, sea level-based N','LineWidth',1,'Color','blue');
load steady_state_coupled.mat
plot(params.sigma_elem*xg/1000,h*params.h0,'DisplayName','Generalized Coulomb Law, coupled hydrology','LineWidth',1,'Color','k','LineStyle',':');
set(ax1,'XTickLabel',[]);
ylabel('Ice thickness [m]','Interpreter','latex')
legend(Location="southwest")

ax2 = subplot('Position',[.1 .62 .8 .17]);
% Plot ice velocity changes with coupling 
load SSA_simple_result.mat
plot(params.sigma*xg*params.xscale/1000,u*params.uscale*params.year,'DisplayName','Power Law','LineWidth',3,'Color','r');
hold on;
load SSA_const_N_result.mat
plot(params.sigma*xg/1000,u*params.u0*params.year,'DisplayName','Generalized Coulomb Law, uniform N','LineWidth',2,'Color','k','LineStyle','--');
load SSA_sea_level_result.mat
plot(params.sigma*xg/1000,u*params.u0*params.year,'DisplayName','Generalized Coulomb Law, sea level-based N','LineWidth',1,'Color','blue');
load steady_state_coupled.mat
plot(params.sigma*xg/1000,u*params.u0*params.year,'DisplayName','Generalized Coulomb Law, coupled hydrology','LineWidth',1,'Color','k','LineStyle',':');
set(ax2,'XTickLabel',[]);
ylabel('Ice velocity [m s$^{-1}$]','Interpreter','latex')

ax3 = subplot('Position',[.1 .43 .8 .17]);
% Plot effective pressure with coupling 
load SSA_const_N_result.mat
plot(params.sigma*xg/1000,100000*ones(size(params.sigma)),'DisplayName','Generalized Coulomb Law, uniform N','LineWidth',2,'Color','k','LineStyle','--');
hold on;
load steady_state_coupled.mat
plot(params.sigma_h*xg/1000,N*params.N0,'DisplayName','Generalized Coulomb Law, coupled hydrology','LineWidth',1,'Color','k','LineStyle',':');
ylabel('Effective pressure [Pa]','Interpreter','latex')
load SSA_sea_level_result.mat
yyaxis("right");
plot(params.sigma*xg/1000,N*params.N0,'DisplayName','Generalized Coulomb Law, sea level-based N','LineWidth',1,'Color','blue');
set(ax3,'XTickLabel',[]);
ax3.YAxis(2).Color = 'blue';
ylabel('Effective pressure [Pa]','Interpreter','latex')


ax4 = subplot('Position',[.1 .24 .8 .17]);
% Plot channel flow with coupling 
load steady_state_coupled.mat
plot(params.sigma_h*xg/1000,Q*params.Q0,'DisplayName','Generalized Coulomb Law, coupled hydrology','LineWidth',1,'Color','k','LineStyle',':');
set(ax4,'XTickLabel',[]);
ylabel('Channel flow [m$^3$ s$^{-1}$]','Interpreter','latex')

ax5 = subplot('Position',[.1 .04 .8 .17]);
% Plot channel flow with coupling 
plot(params.sigma_h*xg/1000,S*params.S0,'DisplayName','Generalized Coulomb Law, coupled hydrology','LineWidth',1,'Color','k','LineStyle',':');
ylabel('Channel area [m$^2$]','Interpreter','latex')
xlabel('Distance from divide [km]','Interpreter','latex')


linkaxes([ax1,ax2,ax3,ax4,ax5],'x');
xlim([150,450]);