load coulomb_steady_state.mat;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
fig = figure(Name='Figure 4',NumberTitle='off');
a1 = subplot('Position',[.1 .8 .8 .15]); 
plot(params.sigma_h,Q,'k','LineWidth',2.5,'DisplayName','S1.C'); ylabel('$Q$','Interpreter','latex');title('Nondimensional steady state variables','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
a2 = subplot('position',[.1 .625 .8 .15]);
plot(params.sigma_h,N,'k','LineWidth',1,'DisplayName','S1.C');ylabel('$N$','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
a3 = subplot('Position',[.1 .45 .8 .15]);
plot(params.sigma_h,S,'k','LineWidth',1,'DisplayName','S1.C');ylabel('$S$','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
a4 = subplot('Position',[.1 .275 .8 .15]);    
plot(params.sigma_elem,h,'k','LineWidth',1,'DisplayName','S1.C');ylabel('$h$','Interpreter','latex');hold on;
set(gca,'Xticklabel',[]);
a5 = subplot('Position',[.1 .1 .8 .15]);
plot(params.sigma,u,'k','LineWidth',1,'DisplayName','S1.C');hold on;ylabel('$u$','Interpreter','latex');xlabel('normalized distance','Interpreter','latex');
linkaxes([a1,a2,a3,a4,a5],'x')

%%
load budd_steady_state.mat;
plot(a1,params.sigma_h,Q, '-r','LineWidth',1,'DisplayName','S1.B');
plot(a2,params.sigma_h,N, '-r','LineWidth',1,'DisplayName','S1.B');
plot(a3,params.sigma_h,S, '-r','LineWidth',1,'DisplayName','S1.B');
plot(a4,params.sigma_elem,h, '-r','LineWidth',1,'DisplayName','S1.B');
plot(a5,params.sigma,u, '-r','LineWidth',1,'DisplayName','S1.B');
%%
load quad_oneway_1.mat;
plot(a1,params.sigma_h,Qs(end,:), '-.k','LineWidth',1,'DisplayName','S2');
plot(a2,params.sigma_h,Ns(end,:), '-.k','LineWidth',1,'DisplayName','S2');
plot(a3,params.sigma_h,Ss(end,:), '-.k','LineWidth',1,'DisplayName','S2');
plot(a4,params.sigma_elem,params.h, '-.k','LineWidth',1,'DisplayName','S2');
plot(a5,params.sigma,params.u, '-.k','LineWidth',1,'DisplayName','S2');

%% S3 numerical solution
%x = linspace(0,1,100);
%N = (1/params.delta).*((x-1) + params.r.*sqrt(1-x));
%plot(a1,x,ones(size(x)),'--r','LineWidth',1,'DisplayName','S3');
%plot(a2,x,N,'-r','LineWidth',1,'DisplayName','S3');
legend(a2,'Orientation','horizontal','Interpreter','latex');
fontsize(fig, 14, "points")

%plot(a3,x,ones(size(x)),'--r','LineWidth',1,'DisplayName','S3');

