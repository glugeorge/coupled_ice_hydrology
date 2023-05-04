clear variables;
%% Load data
load run_1000_full.mat
Q_1000 = Q_dim(end,:);
N_1000 = N_dim(end,:);
S_1000 = S_dim(end,:);
h_1000 = hs(end,:).*params.h0;
u_1000 = us(end,:).*params.u0;
xg_1000 = xgs(end).*params.x0./1000;

load run_100_full.mat
Q_100 = Q_dim(end,:);
N_100 = N_dim(end,:);
S_100 = S_dim(end,:);
h_100 = hs(end,:).*params.h0;
u_100 = us(end,:).*params.u0;
xg_100 = xgs(end).*params.x0./1000;


load run_10_full.mat
Q_10 = Q_dim(end,:);
N_10 = N_dim(end,:);
S_10 = S_dim(end,:);
h_10 = hs(end,:).*params.h0;
u_10 = us(end,:).*params.u0;
xg_10 = xgs(end).*params.x0./1000;

load run_2_full.mat
Q_2= Q_dim(end,:);
N_2= N_dim(end,:);
S_2= S_dim(end,:);
h_2= hs(end,:).*params.h0;
u_2= us(end,:).*params.u0;
xg_2= xgs(end).*params.x0./1000;

%% Plotting against each other

figure("Name","Ice Sheet after 100 years");
title("Grounding line positions");
subplot(3,1,1);
semilogx([2 10 100 1000],[xg_2,xg_10,xg_100,xg_1000],'o');
ylabel('distance (km)');xlabel('Step Size');

subplot(3,1,2)
title("Glacier thicknesses");
plot(params.sigma_elem.*xg_2,h_2,'DisplayName','Nt = 2');hold on;
plot(params.sigma_elem.*xg_10,h_10,'DisplayName','Nt = 10');
plot(params.sigma_elem.*xg_100,h_100,'DisplayName','Nt = 100');
plot(params.sigma_elem.*xg_1000,h_1000,'DisplayName','Nt = 1000');
legend;
xlabel('distance (km)');ylabel('thickness (m)');

subplot(3,1,3);
title("Glacier velocities");
plot(params.sigma.*xg_2,u_2,'DisplayName','Nt = 2');hold on;
plot(params.sigma.*xg_10,u_10,'DisplayName','Nt = 10');
plot(params.sigma.*xg_100,u_100,'DisplayName','Nt = 100');
plot(params.sigma.*xg_1000,u_1000,'DisplayName','Nt = 1000');
legend;
xlabel('distance (km)');ylabel('velocity (m/y)');

figure("Name","Hydrology after 100 years")
subplot(3,1,1);
title("Flow rates");
plot(params.sigma_h.*xg_2,Q_2,'DisplayName','Nt = 2');hold on;
plot(params.sigma_h.*xg_10,Q_10,'DisplayName','Nt = 10');
plot(params.sigma_h.*xg_100,Q_100,'DisplayName','Nt = 100');
plot(params.sigma_h.*xg_1000,Q_1000,'DisplayName','Nt = 1000');
legend;
xlabel('distance (km)');ylabel('flow rate (m^3/s)');

subplot(3,1,2);
title("Effective pressures");
plot(params.sigma_h.*xg_2,N_2,'DisplayName','Nt = 2');hold on;
plot(params.sigma_h.*xg_10,N_10,'DisplayName','Nt = 10');
plot(params.sigma_h.*xg_100,N_100,'DisplayName','Nt = 100');
plot(params.sigma_h.*xg_1000,N_1000,'DisplayName','Nt = 1000');
legend;
xlabel('distance (km)');ylabel('effective pressure (Pa)');

subplot(3,1,3);
title("Surface areas");
plot(params.sigma_h.*xg_2,S_2,'DisplayName','Nt = 2');hold on;
plot(params.sigma_h.*xg_10,S_10,'DisplayName','Nt = 10');
plot(params.sigma_h.*xg_100,S_100,'DisplayName','Nt = 100');
plot(params.sigma_h.*xg_1000,S_1000,'DisplayName','Nt = 1000');
legend;
xlabel('distance (km)');ylabel('area (m^2)');

