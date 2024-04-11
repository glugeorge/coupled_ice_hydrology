%% Loading, analyzing, and plotting transient experiments
clear; close all;
%% Examine Simple-N model to make sure effective pressure profile makes sense
load hybrid_N_coulomb.mat;
params_simpleN = results_hybridN.params;
sigma_simpleN = params_simpleN.sigma;
sigma_elem_simpleN = params_simpleN.sigma_elem;
b = -bed_schoof(results_hybridN.xgs.*sigma_simpleN.*params_simpleN.x0);
h_interp = interp1(sigma_elem_simpleN,results_hybridN.hs.*params_simpleN.h0,sigma_simpleN,'linear','extrap');

N_simpleN = params_simpleN.g.*params_simpleN.rho_i.*h_interp;
% Plotting one slice of this plot(sigma_simpleN,N_simpleN(:,1)) ->
% discontinuity when bed elevation greater than sea level.

%% Plot transient comparison (coulomb)
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
load coulomb_retreat_S.mat;
params = results.params; 
sigma_elem = params.sigma_elem;
sigma_elem = [sigma_elem; sigma_elem(end)];
sigma_h = params.sigma_h;
sigma = params.sigma; 
h = results.hs(:,1);
h = [h; 0];
%figure(Name='Figure 6',NumberTitle='off');

a1 = subplot(2,1,1);
plot(a1,results.ts,results.xgs.*params.x0./1000,'-k','LineWidth',1,'DisplayName','T1.C');hold on; xlabel('model year, $t$ [years]','Interpreter','latex');ylabel('grounding line position, $x_g$ [km]','Interpreter','latex')
ylim([600 1500])
title('Grounding line evolution','Interpreter','latex');
a2 = subplot(2,1,2);
plot(a2,sigma_elem.*results.xgs(1).*params.x0./1000,h.*params.h0+bed_schoof(sigma_elem.*results.xgs(1).*params.x0),'-k','LineWidth',1,'DisplayName','T1.C');ylabel('elevation, h [m]','Interpreter','latex');hold on;
h = results.hs(:,end);
h = [h; 0];
plot(a2,sigma_elem.*results.xgs(end).*params.x0./1000,h.*params.h0+bed_schoof(sigma_elem.*results.xgs(end).*params.x0),'-k','LineWidth',1,'DisplayName','T1.C');ylabel('elevation, h [m]','Interpreter','latex');hold on;
area(a2,linspace(0,1)*1500,bed_schoof(linspace(0,1).*1500e3,params),-2000, 'FaceColor','#808080'); ylim([-2000 6000]);
xlabel('distance from divide, $x$ [km]','Interpreter','latex');
title('Initial and final positions','Interpreter','latex');
%%
load const_N_coulomb.mat;
params = results_constN.params; 
sigma_elem = params.sigma_elem;
sigma_elem = [sigma_elem; sigma_elem(end)];
sigma = params.sigma; 
h = results_constN.hs(:,1);
h = [h; 0];
plot(a1,results_constN.ts,results_constN.xgs.*params.x0./1000,':k','LineWidth',1,'DisplayName','T2.C');
legend(a1,'Orientation','horizontal','Interpreter','latex');

% Plot full connectivity
plot(a1,results_hybridN.ts,results_hybridN.xgs.*params_simpleN.x0./1000,'--k','LineWidth',1,'DisplayName','T3.C')
h = results_hybridN.hs(:,end);
h = [h; 0];
sigma_elem = [sigma_elem_simpleN; sigma_elem_simpleN(end)];
sigma_h = params.sigma_h;
plot(a2,sigma_elem.*results_hybridN.xgs(end).*params_simpleN.x0./1000,h.*params_simpleN.h0+bed_schoof(sigma_elem.*results_hybridN.xgs(end).*params_simpleN.x0),'--k','LineWidth',1,'DisplayName','T1.C');

