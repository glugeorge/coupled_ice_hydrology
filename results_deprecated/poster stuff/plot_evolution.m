set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

figure("Color",'white','Position',[440,136,887,661])

load evolution_to_coupled_hydro_long.mat
ax1 = subplot(3,2,1);
plot(ts,xgs.*params.x0./1e3,'linewidth',3);
xlabel('time [year]','Interpreter','latex');
ylabel('Position [km]',Interpreter='latex');
title('Grounding line',Interpreter='latex');

ax2 = subplot(3,2,2);
contourf(ts,params.sigma_h,Q_dim');colorbar;
xlabel('time (yr)',Interpreter='latex');
ylabel('distance',Interpreter='latex');
title('Flow (m^3/s)',Interpreter='latex');set(gca,'Ydir','Reverse')

ax3 = subplot(3,2,3);contourf(ts,params.sigma_elem,hs'.*params.h0);colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse');
ax4 = subplot(3,2,4);contourf(ts,params.sigma_h,N_dim');colorbar;xlabel('time (yr)');ylabel('distance');title('Effective Pressure (Pa)');set(gca,'Ydir','Reverse')

ax5 = subplot(3,2,5);contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse');
ax6 = subplot(3,2,6);contourf(ts,params.sigma_h,S_dim');colorbar;xlabel('time (yr)');ylabel('distance');title('Surface Area (m^2)');set(gca,'Ydir','Reverse')
