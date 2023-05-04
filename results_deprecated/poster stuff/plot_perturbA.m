set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure("Color",'white','Position',[440,136,887,661])

ax1 = subplot(3,4,1);
load run_2.9_combined.mat
plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time [yr]','Interpreter','latex');ylabel('Grounding line position [km]','Interpreter','latex');title('A = 2.9e-25; Coupled Model','Interpreter','latex')
ax2 = subplot(3,4,5);
contourf(ts,params.sigma_elem,hs'.*params.h0);
c1=colorbar; c1.TickLabelInterpreter="latex"; ylabel(c1,'thickness [m]','Interpreter','latex');
xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
clim([300,2000]);
ax3 = subplot(3,4,9);
contourf(ts,params.sigma,us'.*params.u0.*params.year);
c2 = colorbar; c2.TickLabelInterpreter="latex"; ylabel(c2,'velocity [m yr$^{-1}$]','Interpreter','latex');
xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
%linkaxes([ax1,ax2,ax3],'x');
ax4 = subplot(3,4,2);

load run_2.9_constN.mat
plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time [yr]','Interpreter','latex');ylabel('Grounding line position [km]','Interpreter','latex');title('A = 2.9e-25; Uniform-N Model','Interpreter','latex')
ax5 = subplot(3,4,6);
contourf(ts,params.sigma_elem,hs'.*params.h0);colorbar;xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
c3=colorbar; c3.TickLabelInterpreter="latex"; ylabel(c3,'thickness [m]','Interpreter','latex');
xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
clim([300,2000]);
ax6 = subplot(3,4,10);
contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
c4 = colorbar; c4.TickLabelInterpreter="latex"; ylabel(c4,'velocity [m yr$^{-1}$]','Interpreter','latex');
xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
%linkaxes([ax4,ax5,ax6],'x');

ax7 = subplot(3,4,3);
load run_4.9_combined.mat
plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time [yr]','Interpreter','latex');ylabel('Grounding line position [km]','Interpreter','latex');title('A = 4.9e-25; Coupled Model','Interpreter','latex')
ax8 = subplot(3,4,7);
contourf(ts,params.sigma_elem,hs'.*params.h0);colorbar;xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
c5=colorbar; c5.TickLabelInterpreter="latex"; ylabel(c5,'thickness [m]','Interpreter','latex');
xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
clim([300,2000]);
ax9 = subplot(3,4,11);
contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
c6 = colorbar; c6.TickLabelInterpreter="latex"; ylabel(c6,'velocity [m yr$^{-1}$]','Interpreter','latex');
xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
%linkaxes([ax7,ax8,ax9],'x');

ax10 = subplot(3,4,4);
load run_4.9_constN.mat
plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time [yr]','Interpreter','latex');ylabel('Grounding line position [km]','Interpreter','latex');title('A = 4.9e-25; Uniform-N Model','Interpreter','latex')
ax11 = subplot(3,4,8);
contourf(ts,params.sigma_elem,hs'.*params.h0);colorbar;xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
c7=colorbar; c7.TickLabelInterpreter="latex"; ylabel(c7,'thickness [m]','Interpreter','latex');
xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
clim([300,2000]);
ax12 = subplot(3,4,12);
contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
c8 = colorbar; c8.TickLabelInterpreter="latex"; ylabel(c8,'velocity [m yr$^{-1}$]','Interpreter','latex');
xlabel('time [yr]','Interpreter','latex');ylabel('$\sigma$','Interpreter','latex');set(gca,'Ydir','Reverse');
%linkaxes([ax10,ax11,ax12],'x');
