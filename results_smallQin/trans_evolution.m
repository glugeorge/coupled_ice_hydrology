clear;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
load coulomb_retreat_5000yr.mat
params = results.params; 
sigma_elem = params.sigma_elem;
sigma_elem = [sigma_elem; sigma_elem(end)];
sigma_h = params.sigma_h;
sigma = params.sigma; 
h = results.hs(:,1);
h = [h; 0];
figure(Name='Figure 6',NumberTitle='off');

a1 = subplot(2,1,1);
plot(a1,results.ts,results.xgs.*params.x0./1000,':k','LineWidth',1,'DisplayName','T1.C');hold on; xlabel('model year, $t$ [years]','Interpreter','latex');ylabel('grounding line position, $x_g$ [km]','Interpreter','latex')
ylim([600 1500])
title('Grounding line evolution','Interpreter','latex');
a2 = subplot(2,1,2);
plot(a2,sigma_elem.*results.xgs(1).*params.x0./1000,h.*params.h0+bed_schoof(sigma_elem.*results.xgs(1).*params.x0),':k','LineWidth',1,'DisplayName','T1.C');ylabel('elevation, h [m]','Interpreter','latex');hold on;
h = results.hs(:,end);
h = [h; 0];
plot(a2,sigma_elem.*results.xgs(end).*params.x0./1000,h.*params.h0+bed_schoof(sigma_elem.*results.xgs(end).*params.x0),':k','LineWidth',1,'DisplayName','T1.C');ylabel('elevation, h [m]','Interpreter','latex');hold on;
area(a2,linspace(0,1)*1500,bed_schoof(linspace(0,1).*1500e3,params),-2000, 'FaceColor','#808080'); ylim([-2000 6000]);
xlabel('distance from divide, $x$ [km]','Interpreter','latex');
title('Initial and final steady state positions','Interpreter','latex');
%%
load budd_retreat_5000yr.mat
params = results.params; 
sigma_elem = params.sigma_elem;
sigma_elem = [sigma_elem; sigma_elem(end)];
sigma_h = params.sigma_h;
sigma = params.sigma; 
h = results.hs(:,1);
h = [h; 0];
plot(a1,results.ts,results.xgs.*params.x0./1000,':r','LineWidth',1,'DisplayName','T1.B');
plot(a2,sigma_elem.*results.xgs(1).*params.x0./1000,h.*params.h0+bed_schoof(sigma_elem.*results.xgs(1).*params.x0),':r','LineWidth',1,'DisplayName','T1.B');
h = results.hs(:,end);
h = [h; 0];
plot(a2,sigma_elem.*results.xgs(end).*params.x0./1000,h.*params.h0+bed_schoof(sigma_elem.*results.xgs(end).*params.x0),':r','LineWidth',1,'DisplayName','T1.B');

load const_N_coulomb_5000.mat
params = results_constN.params; 
plot(a1,results.ts,results_constN.xgs.*params.x0./1000,'--k','LineWidth',1,'DisplayName','T2.C');

load const_N_budd_5000.mat
params = results_constN.params; 
plot(a1,results.ts,results_constN.xgs.*params.x0./1000,'--r','LineWidth',1,'DisplayName','T2.B');


legend(a1,'Orientation','horizontal','Interpreter','latex');

