clear variables;
%% Load everything
load run_0_full.mat
Q_0 = Qs.*params.Q0;
N_0 = Ns.*params.N0;
S_0 = Ss.*params.S0;
h_0 = hs.*params.h0;
u_0 = us.*params.u0.*params.year;
xg_0 = xgs.*params.x0./1000;

load run_0.1_full.mat
Q_01 = Qs.*params.Q0;
N_01 = Ns.*params.N0;
S_01 = Ss.*params.S0;
h_01 = hs.*params.h0;
u_01 = us.*params.u0.*params.year;
xg_01 = xgs.*params.x0./1000;

load run_1_full.mat
Q_1 = Qs.*params.Q0;
N_1 = Ns.*params.N0;
S_1 = Ss.*params.S0;
h_1 = hs.*params.h0;
u_1 = us.*params.u0.*params.year;
xg_1 = xgs.*params.x0./1000;

load run_5_full.mat
Q_5 = Qs.*params.Q0;
N_5 = Ns.*params.N0;
S_5 = Ss.*params.S0;
h_5 = hs.*params.h0;
u_5 = us.*params.u0.*params.year;
xg_5 = xgs.*params.x0./1000;

load run_10_full.mat
Q_10 = Qs.*params.Q0;
N_10 = Ns.*params.N0;
S_10 = Ss.*params.S0;
h_10 = hs.*params.h0;
u_10 = us.*params.u0.*params.year;
xg_10 = xgs.*params.x0./1000;

load run_20_full.mat
Q_20 = Qs.*params.Q0;
N_20 = Ns.*params.N0;
S_20 = Ss.*params.S0;
h_20 = hs.*params.h0;
u_20 = us.*params.u0.*params.year;
xg_20 = xgs.*params.x0./1000;

load run_50_full.mat
Q_50 = Qs.*params.Q0;
N_50 = Ns.*params.N0;
S_50 = Ss.*params.S0;
h_50 = hs.*params.h0;
u_50 = us.*params.u0.*params.year;
xg_50 = xgs.*params.x0./1000;

load run_100_full.mat
Q_100 = Qs.*params.Q0;
N_100 = Ns.*params.N0;
S_100 = Ss.*params.S0;
h_100 = hs.*params.h0;
u_100 = us.*params.u0.*params.year;
xg_100 = xgs.*params.x0./1000;

load run_500_full.mat
Q_500 = Qs.*params.Q0;
N_500 = Ns.*params.N0;
S_500 = Ss.*params.S0;
h_500 = hs.*params.h0;
u_500 = us.*params.u0.*params.year;
xg_500 = xgs.*params.x0./1000;



%% Plot everything after 500 years 
figure("Name","Ice Sheet after 500 years");
title("Grounding line positions");
subplot(3,1,1);
semilogx([0 0.1 1 5 10 20 50 100 500],[xg_0(end) xg_01(end) xg_1(end) xg_5(end) xg_10(end) xg_20(end) xg_50(end) xg_100(end) xg_500(end)],'o');
xlabel('Flow rate (m^3/s)');ylabel('Grounding line position (km)')

subplot(3,1,2)
title("Glacier thicknesses");
plot(params.sigma_elem.*xg_0(end),h_0(end,:),'DisplayName','Q = 0');hold on;
plot(params.sigma_elem.*xg_01(end),h_01(end,:),'DisplayName','Q = 0.1');
plot(params.sigma_elem.*xg_1(end),h_1(end,:),'DisplayName','Q = 1');
plot(params.sigma_elem.*xg_5(end),h_5(end,:),'DisplayName','Q = 5');
plot(params.sigma_elem.*xg_10(end),h_10(end,:),'DisplayName','Q = 10');
plot(params.sigma_elem.*xg_20(end),h_20(end,:),'DisplayName','Q = 20');
plot(params.sigma_elem.*xg_50(end),h_50(end,:),'DisplayName','Q = 50');
plot(params.sigma_elem.*xg_100(end),h_100(end,:),'DisplayName','Q = 100');
plot(params.sigma_elem.*xg_500(end),h_500(end,:),'DisplayName','Q = 500');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("thickness (m)");

subplot(3,1,3);
title("Glacier velocities");
plot(params.sigma_elem.*xg_0(end),u_0(end,:),'DisplayName','Q = 0');hold on;
plot(params.sigma_elem.*xg_01(end),u_01(end,:),'DisplayName','Q = 0.1');
plot(params.sigma_elem.*xg_1(end),u_1(end,:),'DisplayName','Q = 1');
plot(params.sigma_elem.*xg_5(end),u_5(end,:),'DisplayName','Q = 5');
plot(params.sigma_elem.*xg_10(end),u_10(end,:),'DisplayName','Q = 10');
plot(params.sigma_elem.*xg_20(end),u_20(end,:),'DisplayName','Q = 20');
plot(params.sigma_elem.*xg_50(end),u_50(end,:),'DisplayName','Q = 50');
plot(params.sigma_elem.*xg_100(end),u_100(end,:),'DisplayName','Q = 100');
plot(params.sigma_elem.*xg_500(end),u_500(end,:),'DisplayName','Q = 500');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("velocity (m/y)");

figure("Name","Hydrology after 500 years")
subplot(3,1,1);
title("Flow rates");
plot(params.sigma_elem.*xg_0(end),Q_0(end,:),'DisplayName','Q = 0');hold on;
plot(params.sigma_elem.*xg_01(end),Q_01(end,:),'DisplayName','Q = 0.1');
plot(params.sigma_elem.*xg_1(end),Q_1(end,:),'DisplayName','Q = 1');
plot(params.sigma_elem.*xg_5(end),Q_5(end,:),'DisplayName','Q = 5');
plot(params.sigma_elem.*xg_10(end),Q_10(end,:),'DisplayName','Q = 10');
plot(params.sigma_elem.*xg_20(end),Q_20(end,:),'DisplayName','Q = 20');
plot(params.sigma_elem.*xg_50(end),Q_50(end,:),'DisplayName','Q = 50');
plot(params.sigma_elem.*xg_100(end),Q_100(end,:),'DisplayName','Q = 100');
plot(params.sigma_elem.*xg_500(end),Q_500(end,:),'DisplayName','Q = 500');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("flow (m^3/s)");


subplot(3,1,2);
title("Effective pressures");
plot(params.sigma_elem.*xg_0(end),N_0(end,:),'DisplayName','Q = 0');hold on;
plot(params.sigma_elem.*xg_01(end),N_01(end,:),'DisplayName','Q = 0.1');
plot(params.sigma_elem.*xg_1(end),N_1(end,:),'DisplayName','Q = 1');
plot(params.sigma_elem.*xg_5(end),N_5(end,:),'DisplayName','Q = 5');
plot(params.sigma_elem.*xg_10(end),N_10(end,:),'DisplayName','Q = 10');
plot(params.sigma_elem.*xg_20(end),N_20(end,:),'DisplayName','Q = 20');
plot(params.sigma_elem.*xg_50(end),N_50(end,:),'DisplayName','Q = 50');
plot(params.sigma_elem.*xg_100(end),N_100(end,:),'DisplayName','Q = 100');
plot(params.sigma_elem.*xg_500(end),N_500(end,:),'DisplayName','Q = 500');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("effective pressure (Pa)");

subplot(3,1,3);
title("Surface areas");
plot(params.sigma_elem.*xg_0(end),S_0(end,:),'DisplayName','Q = 0');hold on;
plot(params.sigma_elem.*xg_01(end),S_01(end,:),'DisplayName','Q = 0.1');
plot(params.sigma_elem.*xg_1(end),S_1(end,:),'DisplayName','Q = 1');
plot(params.sigma_elem.*xg_5(end),S_5(end,:),'DisplayName','Q = 5');
plot(params.sigma_elem.*xg_10(end),S_10(end,:),'DisplayName','Q = 10');
plot(params.sigma_elem.*xg_20(end),S_20(end,:),'DisplayName','Q = 20');
plot(params.sigma_elem.*xg_50(end),S_50(end,:),'DisplayName','Q = 50');
plot(params.sigma_elem.*xg_100(end),S_100(end,:),'DisplayName','Q = 100');
plot(params.sigma_elem.*xg_500(end),S_500(end,:),'DisplayName','Q = 500');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("area (m^2)");

