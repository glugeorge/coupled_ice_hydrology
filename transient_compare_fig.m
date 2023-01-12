close all; clear;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
load schoof_retreat_highres_tiny.mat
load const_N.mat
results_constcoulomb = results_constN;
load const_powerlaw.mat;
results_constpower = results_constN;

tiledlayout(4,1, 'TileSpacing', 'tight'); 

nexttile;
plot(results_constcoulomb.ts,results_constcoulomb.xgs.*results_constcoulomb.params.x0./1e3,'linewidth',1,'DisplayName','Coulomb Law');hold on;
plot(results_constpower.ts,results_constpower.xgs.*results_constpower.params.x0./1e3,'--','linewidth',1,'DisplayName','Power Law');
xlabel('time [yr]','Interpreter','latex');ylabel('$x_g$ [km]','Interpreter','latex');
plot(results.ts,results.xgs.*results.params.x0./1e3,'linewidth',1,'DisplayName','Coupled')

nexttile;
plot_profile(results_constcoulomb,'Coulomb Law',results_constcoulomb.hs(:,end).*results_constcoulomb.params.h0,'h [m]','-',results_constcoulomb.params.sigma_elem);hold on;
plot_profile(results_constpower,'Power Law',results_constpower.hs(:,end).*results_constpower.params.h0,'h [m]','--',results_constpower.params.sigma_elem)
plot_profile(results,'Coupled',results.hs(:,end).*results.params.h0,'h [m]','-',results.params.sigma_elem);
xlim([1300 1375]); set(gca,'Xticklabel',[]);

nexttile;
plot_profile(results_constcoulomb,'Coulomb Law',results_constcoulomb.us(:,end).*results_constcoulomb.params.u0*60*60*24*365,'u [m y$^{-1}$]','-',results_constcoulomb.params.sigma);hold on;
plot_profile(results_constpower,'Power Law',results_constpower.us(:,end).*results_constpower.params.u0*60*60*24*365,'u [m y$^{-1}$]','--',results_constpower.params.sigma)
plot_profile(results,'Coupled',results.us(:,end).*results.params.u0*60*60*24*365,'u [m y$^{-1}$]','-',results.params.sigma);
xlim([1300 1375]); set(gca,'Xticklabel',[]);legend('Orientation','horizontal','Interpreter','latex');

nexttile;
plot(results.params.sigma_h.*results.xgs(1).*results.params.x0./1000,results.Ns(:,1).*results.params.N0,'-','LineWidth',1,'DisplayName','Static N');hold on;
plot(results.params.sigma_h.*results.xgs(1).*results.params.x0./1000,results.Ns(:,1).*results.params.N0,'--','LineWidth',1,'DisplayName','Static N');
plot_profile(results,'Coupled N',results.Ns(:,end).*results.params.N0,'N [Pa]','-',results.params.sigma_h);
xlim([1300 1375]);xlabel('distance from divide, $x$ [km]','Interpreter','latex');

function plot_profile(result,result_name,variable,variable_name,linestyle,grid);
    params = result.params;
    grid_dim = grid.*result.xgs(end).*params.x0./1000;
    plot(grid_dim,variable,linestyle,'LineWidth',1,'DisplayName',result_name);ylabel(variable_name,'Interpreter','latex')

end
