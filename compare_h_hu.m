% Evaluating impact of a realistic h and u imposed on a hydrology model
clear variables
load realistic_h.mat
load realistic_hu.mat
load steady_params.mat
load steady_xg.mat

h_grid = params.sigma_elem*xg;
x_grid = params.x0*linspace(0,h_grid(end),params.Nx)'./1000;% get length of glacier in km


N_diff = N_dim-N_dim_1;
Q_diff = Q_dim-Q_dim_1;
S_diff = S_dim-S_dim_1;

subplot(3,1,1);surface(ts,x_grid,Q_diff',EdgeColor='None');caxis([-1 1]);colorbar;xlabel('time (yr)');ylabel('distance');title('Flow (m^3/s)');set(gca,'Ydir','Reverse')
subplot(3,1,2);surface(ts,x_grid,N_diff',EdgeColor='None');caxis([-4e5 4e5]);colorbar;xlabel('time (yr)');ylabel('distance');title('Effective Pressure (Pa)');set(gca,'Ydir','Reverse')
subplot(3,1,3);surface(ts,x_grid,S_diff',EdgeColor='None');caxis([-5 5]);colorbar;xlabel('time (yr)');ylabel('distance');title('Surface Area (m^2)');set(gca,'Ydir','Reverse')
colormap(bluewhitered);
figure('Name','Evolution of Drainage');
subplot(3,1,1);
plot(ts,N_diff(:,1),'DisplayName','Channel Entrance');
hold on;
plot(ts,N_diff(:,end),'--','DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)','Interpreter','latex');
ylabel('Effective Pressure, \emph{N} (Pa)','Interpreter','latex');
title('Effective Pressure over time');

subplot(3,1,2);
plot(ts,Q_diff(:,1),'DisplayName','Channel Entrance');
hold on;
plot(ts,Q_diff(:,end),'DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)','Interpreter','latex');
ylabel('channel flow rate, \emph{Q} $(m^{3} s^{-1})$','Interpreter','latex');
title('Flow over time');
legend('Location','northwest');

subplot(3,1,3);
plot(ts,S_diff(:,1),'-o','DisplayName','Channel Entrance');
hold on;
plot(ts,S_diff(:,end),'-o','DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)','Interpreter','latex');
ylabel('channel cross-sectional area, \emph{S} $(m^2)$','Interpreter','latex');
title('Channel area over time');
legend('Location','northwest');

figure()
plot(x_grid,u_interp.*params.u0*params.year);
xlabel('distance, \emph{x} (km)','Interpreter','latex');ylabel('velocity, \emph{u} ($m y^{-1}$)','Interpreter','latex');
