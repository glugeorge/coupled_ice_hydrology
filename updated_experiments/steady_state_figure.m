load S1C_highres.mat;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
fig = figure(Name='Figure 4',NumberTitle='off');
a1 = subplot('Position',[.1 .8 .8 .15]); 
plot(params.sigma_h,Q*params.Q0,'k','LineWidth',1,'DisplayName','S1.C'); ylabel('$Q$ [m$^3$s$^{-1}$]','Interpreter','latex');title('Nondimensional steady state variables','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
a2 = subplot('position',[.1 .625 .8 .15]);
plot(params.sigma_h,N*params.N0./1e6,'k','LineWidth',1,'DisplayName','S1.C');ylabel('$N$ [MPa]','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
a3 = subplot('Position',[.1 .45 .8 .15]);
plot(params.sigma_h,S*params.S0,'k','LineWidth',1,'DisplayName','S1.C');ylabel('$S$ [m$^2$]','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
a4 = subplot('Position',[.1 .275 .8 .15]);    
plot(params.sigma_elem,h*params.h0,'k','LineWidth',1,'DisplayName','S1.C');ylabel('$h$ [m]','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
a5 = subplot('Position',[.1 .1 .8 .15]);
plot(params.sigma,u*params.u0*params.year,'k','LineWidth',1,'DisplayName','S1.C');hold on;ylabel('$u$ [m y$^{-1}$]','Interpreter','latex');xlabel('normalized distance','Interpreter','latex');
linkaxes([a1,a2,a3,a4,a5],'x')

%%
load S1B_highres.mat;
plot(a1,params.sigma_h,Q*params.Q0, '-r','LineWidth',1,'DisplayName','S1.B');
plot(a2,params.sigma_h,N*params.N0./1e6, '-r','LineWidth',1,'DisplayName','S1.B');
plot(a3,params.sigma_h,S*params.S0, '-r','LineWidth',1,'DisplayName','S1.B');
plot(a4,params.sigma_elem,h*params.h0, '-r','LineWidth',1,'DisplayName','S1.B');
plot(a5,params.sigma,u*params.u0*params.year, '-r','LineWidth',1,'DisplayName','S1.B');
%%
load S2.mat;
plot(a1,params.sigma_h,Qs(end,:)*params.Q0, '-.k','LineWidth',1,'DisplayName','S2');
plot(a2,params.sigma_h,Ns(end,:)*params.N0./1e6, '-.k','LineWidth',1,'DisplayName','S2');
plot(a3,params.sigma_h,Ss(end,:)*params.S0, '-.k','LineWidth',1,'DisplayName','S2');
plot(a4,params.sigma_elem,params.h*params.h0, '-.k','LineWidth',1,'DisplayName','S2');
plot(a5,params.sigma,params.u*params.u0*params.year, '-.k','LineWidth',1,'DisplayName','S2');

%% S3 numerical solution
%x = linspace(0,1,100);
%N = (1/params.delta).*((x-1) + params.r.*sqrt(1-x));
%plot(a1,x,ones(size(x)),'--r','LineWidth',1,'DisplayName','S3');
%plot(a2,x,N,'-r','LineWidth',1,'DisplayName','S3');
legend(a2,'Orientation','horizontal','Interpreter','latex');
fontsize(fig, 14, "points")

%plot(a3,x,ones(size(x)),'--r','LineWidth',1,'DisplayName','S3');

