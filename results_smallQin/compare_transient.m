clear;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
load coulomb_retreat_50yr.mat
results_c = results;
load budd_retreat_50yr.mat
results_b = results;
load const_N_coulomb_50.mat
results_constc = results_constN;
load const_N_budd_50.mat;
results_constb = results_constN;
figure(Name='Figure 7',NumberTitle='off')
t = tiledlayout(4,2, 'TileSpacing', 'tight'); 
nexttile;
plot(results_c.ts(1:200),results_c.xgs(1:200).*results_c.params.x0./1e3,'-k','linewidth',1,'DisplayName','T1.C');hold on;
plot(results_constc.ts(1:200),results_constc.xgs(1:200).*results_constc.params.x0./1e3,'--k','linewidth',1,'DisplayName','T2.C')
xlabel('time (yr)','Interpreter','latex');ylabel('grounding line position, $x_g$ (km)','Interpreter','latex');
title('Regularized Coulomb law','Interpreter','latex')
legend('Orientation','horizontal','Interpreter','latex');

nexttile;
plot(results_b.ts(1:200),results_b.xgs(1:200).*results_b.params.x0./1e3,'-r','linewidth',1,'DisplayName','T1.B');hold on;
xlabel('time (yr)','Interpreter','latex');
plot(results_constb.ts(1:200),results_constb.xgs(1:200).*results_constb.params.x0./1e3,'--r','linewidth',1,'DisplayName','T2.B');
legend('Orientation','horizontal','Interpreter','latex');
title('Budd law','Interpreter','latex')

nexttile;
plot_profile(results_constc,'T2.C',results_constc.hs(:,200).*results_constc.params.h0,'h (m)','--k',results_constc.params.sigma_elem);hold on;
plot_profile(results_c,'T1.C',results_c.hs(:,200).*results_c.params.h0,'ice thickness, h (m)','-k',results_c.params.sigma_elem);
xline(results_constc.xgs(200).*results_c.params.x0./1e3,':','linewidth',1,Color='black')
xlim([1250 1390]);set(gca,'Xticklabel',[]);

nexttile;
plot_profile(results_constb,'T2.B',results_constb.hs(:,200).*results_constb.params.h0,'h [m]','--r',results_constb.params.sigma_elem); hold on;
plot_profile(results_b,'T1.B',results_b.hs(:,200).*results_b.params.h0,'','-r',results_b.params.sigma_elem);
xline(results_constb.xgs(200).*results_b.params.x0./1e3,':','linewidth',1,Color='black')
xlim([1230 1300]); set(gca,'Xticklabel',[]);


nexttile;
plot_profile(results_c,'T1.C',results_c.us(:,200).*results_c.params.u0*60*60*24*365,'ice velocity, u (m y$^{-1}$)','-k',results_c.params.sigma);hold on;
plot_profile(results_constc,'T2.C',results_constc.us(:,200).*results_constc.params.u0*60*60*24*365,'ice velocity, u (m y$^{-1}$)','--k',results_constc.params.sigma);
xline(results_constc.xgs(200).*results_c.params.x0./1e3,':','linewidth',1,Color='black')
xlim([1250 1390]); set(gca,'Xticklabel',[]);

nexttile;
plot_profile(results_b,'T1.B',results_b.us(:,200).*results_b.params.u0*60*60*24*365,'','-r',results_b.params.sigma);hold on;
plot_profile(results_constb,'T2.B',results_constb.us(:,200).*results_constb.params.u0*60*60*24*365,'','--r',results_constb.params.sigma)
xline(results_constb.xgs(200).*results_b.params.x0./1e3,':','linewidth',1,Color='black')
xlim([1230 1300]); set(gca,'Xticklabel',[]);

nexttile;
plot(results_c.params.sigma_h.*results_c.xgs(1).*results_c.params.x0./1000,results_c.Ns(:,1).*results_c.params.N0,'--k','LineWidth',1,'DisplayName','T2.C');hold on;
plot_profile(results_c,'T1.C',results_c.Ns(:,200).*results_c.params.N0,'effective pressure, N (Pa)','-k',results_c.params.sigma_h);
xline(results_constc.xgs(200).*results_c.params.x0./1e3,':','linewidth',1,Color='black')
N_c_interp = interp1(results_c.params.sigma_h.*results_c.xgs(1).*results_c.params.x0./1000,results_c.Ns(:,1).*results_c.params.N0,results_constc.xgs(200).*results_c.params.x0./1e3);
plot(results_constc.xgs(200).*results_c.params.x0./1e3,N_c_interp,'kx',MarkerSize=10);
xlim([1250 1390]);

nexttile;
plot(results_b.params.sigma_h.*results_b.xgs(1).*results_b.params.x0./1000,results_b.Ns(:,1).*results_b.params.N0,'--r','LineWidth',1,'DisplayName','T2.B');hold on;
plot_profile(results_b,'T1.B',results_b.Ns(:,200).*results_b.params.N0,'','-r',results_b.params.sigma_h);
xline(results_constb.xgs(200).*results_b.params.x0./1e3,':','linewidth',1,Color='black')
N_b_interp = interp1(results_b.params.sigma_h.*results_b.xgs(1).*results_b.params.x0./1000,results_b.Ns(:,1).*results_b.params.N0,results_constb.xgs(200).*results_b.params.x0./1e3);
plot(results_constb.xgs(200).*results_b.params.x0./1e3,N_b_interp,'rx',MarkerSize=10);
xlim([1230 1300]);


xlabel(t,'distance from divide, $x$ (km)','Interpreter','latex');


function plot_profile(result,result_name,variable,variable_name,linestyle,grid)
    params = result.params;
    grid_dim = grid.*result.xgs(200).*params.x0./1000;
    plot(grid_dim,variable,linestyle,'LineWidth',1,'DisplayName',result_name);ylabel(variable_name,'Interpreter','latex')

end
