load fig1_steady_state_combined.mat;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
figure()
ax1 = subplot('Position',[.1 .8 .8 .15]); 
plot(params.sigma_h,Q,'k','LineWidth',1,'DisplayName','Coupled'); ylabel('Q','Interpreter','latex');title('Nondimensionalized steady state variables','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
ax2 = subplot('position',[.1 .625 .8 .15]);
plot(params.sigma_h,N,'k','LineWidth',1,'DisplayName','Coupled');ylabel('N','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
ax3 = subplot('Position',[.1 .45 .8 .15]);
plot(params.sigma_h,S,'k','LineWidth',1,'DisplayName','Coupled');ylabel('S','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
ax4 = subplot('Position',[.1 .275 .8 .15]);    
plot(params.sigma_elem,h,'k','LineWidth',1,'DisplayName','Coupled');ylabel('h','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
ax5 = subplot('Position',[.1 .1 .8 .15]);
plot(params.sigma,u,'k','LineWidth',1,'DisplayName','Coupled');hold on;ylabel('u','Interpreter','latex');xlabel('Normalized distance to grounding line','Interpreter','latex');ylim([0 45])
linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
%%
load slab_oneway.mat;
plot(ax1,params.sigma_h,Qs(end,:), '--b','LineWidth',1,'DisplayName','Slab');
plot(ax2,params.sigma_h,Ns(end,:), '--b','LineWidth',1,'DisplayName','Slab');
plot(ax3,params.sigma_h,Ss(end,:), '--b','LineWidth',1,'DisplayName','Slab');
plot(ax4,params.sigma_elem,params.h, '--b','LineWidth',1,'DisplayName','Slab');
plot(ax5,params.sigma,params.u, '--b','LineWidth',1,'DisplayName','Slab');
%%
load quad_oneway.mat;
plot(ax1,params.sigma_h,Qs(end,:), ':r','LineWidth',1,'DisplayName','Quadratic');
plot(ax2,params.sigma_h,Ns(end,:), ':r','LineWidth',1,'DisplayName','Quadratic');
plot(ax3,params.sigma_h,Ss(end,:), ':r','LineWidth',1,'DisplayName','Quadratic');
plot(ax4,params.sigma_elem,params.h, ':r','LineWidth',1,'DisplayName','Quadratic');
plot(ax5,params.sigma,params.u, ':r','LineWidth',1,'DisplayName','Quadratic');
legend('Orientation','horizontal','Interpreter','latex');

