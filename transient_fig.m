
close all; clear;
set(groot,'defaultAxesTickLabelInterpreter','latex');  

load schoof_retreat_highres_long.mat
params = results.params; 
sigma_elem = params.sigma_elem;
sigma_elem = [sigma_elem; sigma_elem(end)];
sigma_h = params.sigma_h;
sigma = params.sigma; 
tiledlayout(2,1, 'TileSpacing', 'compact'); 
nexttile;
for t=1:10:100
    Q = results.Qs(:,t);
    N = results.Ns(:,t);
    S = results.Ss(:,t);
    h = results.hs(:,t);
    h = [h; 0];
    u = results.us(:,t);
    xg = results.xgs(t);
    if t ==1
        plot(sigma_elem.*xg.*params.x0./1000,h*params.h0+bed_schoof(sigma_elem.*xg.*params.x0),'-k','LineWidth',2);ylabel('elevation, h [m]','Interpreter','latex');hold on;
    else
        plot(sigma_elem.*xg.*params.x0./1000,h*params.h0+bed_schoof(sigma_elem.*xg.*params.x0),'-k');ylabel('elevation, h [m]','Interpreter','latex');hold on;
    end
end
plot(linspace(0,1)*1500,bed_schoof(linspace(0,1).*1500e3,params),'-k'); set(gca,'Xticklabel',[]); ylim([-2000 6000]);


load schoof_advance_highres_long.mat
params = results.params; 
sigma_elem = params.sigma_elem;
sigma_elem = [sigma_elem; sigma_elem(end)];
sigma_h = params.sigma_h;
sigma = params.sigma; 

for t=1:10:100
    Q = results.Qs(:,t);
    N = results.Ns(:,t);
    S = results.Ss(:,t);
    h = results.hs(:,t);
    h = [h; 0];
    u = results.us(:,t);
    xg = results.xgs(t);
    if t ==1
        plot(sigma_elem.*xg.*params.x0./1000,h*params.h0+bed_schoof(sigma_elem.*xg.*params.x0),'--k','LineWidth',2);ylabel('elevation, h [m]','Interpreter','latex');hold on;
        nexttile;
        plot(sigma_elem.*xg.*params.x0./1000,h*params.h0+bed_schoof(sigma_elem.*xg.*params.x0),'--k','LineWidth',2);ylabel('elevation, h [m]','Interpreter','latex');hold on;
    else
        
        plot(sigma_elem.*xg.*params.x0./1000,h*params.h0+bed_schoof(sigma_elem.*xg.*params.x0),'-k');ylabel('elevation, h [m]','Interpreter','latex');hold on;
    end
end
plot(linspace(0,1)*1500,bed_schoof(linspace(0,1).*1500e3,params),'-k'); ylim([-2000 6000]);
xlabel('distance from divide, x [km]','Interpreter','latex');


