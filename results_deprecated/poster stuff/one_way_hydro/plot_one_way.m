clear all;
load realistic_h.mat
load realistic_hu.mat
load SSA_simple_result.mat

h_grid = params.sigma_elem*xg;
x_grid = params.xscale*linspace(0,h_grid(end),params.Nx)'./1000;% get length of glacier in km

figure('Name','Evolution of Drainage');
t = tiledlayout(4,2);

ax1 = nexttile;
contourf(ts,x_grid,Q_dim');
c1 = colorbar; ylabel(c1,'Water Flow [m$^3$ s$^{-1}$]','Interpreter','latex');c1.TickLabelInterpreter="latex";
xlabel('Time [yr]','Interpreter','latex');ylabel('Distance from divide [km]','Interpreter','latex');set(gca,'Ydir','Reverse')
title('One way hydrology with realistic glacier height','Interpreter','latex');
ax2 = nexttile;
contourf(ts,x_grid,Q_dim_1');
c2 = colorbar; ylabel(c2,'Water Flow [m$^3$ s$^{-1}$]','Interpreter','latex');c2.TickLabelInterpreter="latex";
xlabel('Time [yr]','Interpreter','latex');ylabel('Distance from divide [km]','Interpreter','latex');set(gca,'Ydir','Reverse')
title('One way hydrology with realistic glacier height and velocity','Interpreter','latex');

ax3 = nexttile;
contourf(ts,x_grid,N_dim');
c3 = colorbar; ylabel(c3,'Effective Pressure [Pa]','Interpreter','latex');c3.TickLabelInterpreter="latex";
xlabel('Time [yr]','Interpreter','latex');ylabel('Distance from divide [km]','Interpreter','latex');set(gca,'Ydir','Reverse')

ax4 = nexttile;
contourf(ts,x_grid,N_dim_1');
c4 = colorbar; ylabel(c4,'Effective Pressure [Pa]','Interpreter','latex');c4.TickLabelInterpreter="latex";
xlabel('Time [yr]','Interpreter','latex');ylabel('Distance from divide [km]','Interpreter','latex');set(gca,'Ydir','Reverse')

ax5 = nexttile;
surface(ts,x_grid,S_dim',EdgeColor='None');
c5 = colorbar; ylabel(c5,'Channel Area [m$^2$]','Interpreter','latex');c5.TickLabelInterpreter="latex";
xlabel('Time [yr]','Interpreter','latex');ylabel('Distance from divide [km]','Interpreter','latex');set(gca,'Ydir','Reverse')
clim([5,10]);
xlim([ts(1) ts(end)])
ylim([x_grid(1) x_grid(end)])

ax6 = nexttile;
surface(ts,x_grid,S_dim_1',EdgeColor='None');
c5 = colorbar; ylabel(c5,'Channel Area [m$^2$]','Interpreter','latex');c6.TickLabelInterpreter="latex";
xlabel('Time [yr]','Interpreter','latex');ylabel('Distance from divide [km]','Interpreter','latex');set(gca,'Ydir','Reverse')
clim([5,10]);
xlim([ts(1) ts(end)])
ylim([x_grid(1) x_grid(end)])

nexttile([1,2])
plot(ts,S_dim(:,end),'-o','DisplayName','No ice velocity');
hold on;
plot(ts,S_dim_1(:,end),'-o','DisplayName','With ice velocity');
xlabel('Time [yr]','Interpreter','latex');
ylabel('Channel Area [m$^2$]','Interpreter','latex');
title('Evolution of channel area at terminus','Interpreter','latex');
legend('Location','northwest');

linkaxes([ax1,ax3,ax5],'x')