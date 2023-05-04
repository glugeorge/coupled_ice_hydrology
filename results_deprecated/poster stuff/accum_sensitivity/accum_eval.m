clear variables;
%% Load everything
load run_0.0_full.mat
Q_0 = Qs.*params.Q0;
N_0 = Ns.*params.N0;
S_0 = Ss.*params.S0;
h_0 = hs.*params.h0;
u_0 = us.*params.u0.*params.year;
xg_0 = xgs.*params.x0./1000;

load run_0.2_full.mat
Q_02 = Qs.*params.Q0;
N_02 = Ns.*params.N0;
S_02 = Ss.*params.S0;
h_02 = hs.*params.h0;
u_02 = us.*params.u0.*params.year;
xg_02 = xgs.*params.x0./1000;

load run_0.4_full.mat
Q_04 = Qs.*params.Q0;
N_04 = Ns.*params.N0;
S_04 = Ss.*params.S0;
h_04 = hs.*params.h0;
u_04 = us.*params.u0.*params.year;
xg_04 = xgs.*params.x0./1000;

load run_0.6_full.mat
Q_06 = Qs.*params.Q0;
N_06 = Ns.*params.N0;
S_06 = Ss.*params.S0;
h_06 = hs.*params.h0;
u_06 = us.*params.u0.*params.year;
xg_06 = xgs.*params.x0./1000;

load run_0.8_full.mat
Q_08 = Qs.*params.Q0;
N_08 = Ns.*params.N0;
S_08 = Ss.*params.S0;
h_08 = hs.*params.h0;
u_08 = us.*params.u0.*params.year;
xg_08 = xgs.*params.x0./1000;

load run_1.0_full.mat
Q_1 = Qs.*params.Q0;
N_1 = Ns.*params.N0;
S_1 = Ss.*params.S0;
h_1 = hs.*params.h0;
u_1 = us.*params.u0.*params.year;
xg_1 = xgs.*params.x0./1000;

load run_1.2_full.mat
Q_12 = Qs.*params.Q0;
N_12 = Ns.*params.N0;
S_12 = Ss.*params.S0;
h_12 = hs.*params.h0;
u_12 = us.*params.u0.*params.year;
xg_12 = xgs.*params.x0./1000;

load run_1.4_full.mat
Q_14 = Qs.*params.Q0;
N_14 = Ns.*params.N0;
S_14 = Ss.*params.S0;
h_14 = hs.*params.h0;
u_14 = us.*params.u0.*params.year;
xg_14 = xgs.*params.x0./1000;

load run_1.6_full.mat
Q_16 = Qs.*params.Q0;
N_16 = Ns.*params.N0;
S_16 = Ss.*params.S0;
h_16 = hs.*params.h0;
u_16 = us.*params.u0.*params.year;
xg_16 = xgs.*params.x0./1000;

load run_1.8_full.mat
Q_18 = Qs.*params.Q0;
N_18 = Ns.*params.N0;
S_18 = Ss.*params.S0;
h_18 = hs.*params.h0;
u_18 = us.*params.u0.*params.year;
xg_18 = xgs.*params.x0./1000;

load run_2.0_full.mat
Q_2 = Qs.*params.Q0;
N_2 = Ns.*params.N0;
S_2 = Ss.*params.S0;
h_2 = hs.*params.h0;
u_2 = us.*params.u0.*params.year;
xg_2 = xgs.*params.x0./1000;

%% Plot everything after 500 years 
figure("Name","Ice Sheet after 500 years");
title("Grounding line positions");
subplot(3,1,1);
plot([0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2],[xg_0(end) xg_02(end) xg_04(end) xg_06(end) xg_08(end) xg_1(end) xg_12(end) xg_14(end) xg_16(end) xg_18(end) xg_2(end)],'o');
xlabel('accumulation rate (m/y)');ylabel("grounding line posiiton (km)");

subplot(3,1,2)
title("Glacier thicknesses");
plot(params.sigma_elem.*xg_0(end),h_0(end,:),'DisplayName','a = 0');hold on;
plot(params.sigma_elem.*xg_02(end),h_02(end,:),'DisplayName','a = 0.2');
plot(params.sigma_elem.*xg_04(end),h_04(end,:),'DisplayName','a = 0.4');
plot(params.sigma_elem.*xg_06(end),h_06(end,:),'DisplayName','a = 0.6');
plot(params.sigma_elem.*xg_08(end),h_08(end,:),'DisplayName','a = 0.8');
plot(params.sigma_elem.*xg_1(end),h_1(end,:),'DisplayName','a = 1');
plot(params.sigma_elem.*xg_12(end),h_12(end,:),'DisplayName','a = 1.2');
plot(params.sigma_elem.*xg_14(end),h_14(end,:),'DisplayName','a = 1.4');
plot(params.sigma_elem.*xg_16(end),h_16(end,:),'DisplayName','a = 1.6');
plot(params.sigma_elem.*xg_18(end),h_18(end,:),'DisplayName','a = 1.8');
plot(params.sigma_elem.*xg_2(end),h_2(end,:),'DisplayName','a = 2');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("thickness (m)");

subplot(3,1,3);
title("Glacier velocities");
plot(params.sigma_elem.*xg_0(end),u_0(end,:),'DisplayName','a = 0');hold on;
plot(params.sigma_elem.*xg_02(end),u_02(end,:),'DisplayName','a = 0.2');
plot(params.sigma_elem.*xg_04(end),u_04(end,:),'DisplayName','a = 0.4');
plot(params.sigma_elem.*xg_06(end),u_06(end,:),'DisplayName','a = 0.6');
plot(params.sigma_elem.*xg_08(end),u_08(end,:),'DisplayName','a = 0.8');
plot(params.sigma_elem.*xg_1(end),u_1(end,:),'DisplayName','a = 1');
plot(params.sigma_elem.*xg_12(end),u_12(end,:),'DisplayName','a = 1.2');
plot(params.sigma_elem.*xg_14(end),u_14(end,:),'DisplayName','a = 1.4');
plot(params.sigma_elem.*xg_16(end),u_16(end,:),'DisplayName','a = 1.6');
plot(params.sigma_elem.*xg_18(end),u_18(end,:),'DisplayName','a = 1.8');
plot(params.sigma_elem.*xg_2(end),u_2(end,:),'DisplayName','a = 2');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("velocity (m/y)");

figure("Name","Hydrology after 500 years")
subplot(3,1,1);
title("Flow rates");
plot(params.sigma_elem.*xg_0(end),Q_0(end,:),'DisplayName','a = 0');hold on;
plot(params.sigma_elem.*xg_02(end),Q_02(end,:),'DisplayName','a = 0.2');
plot(params.sigma_elem.*xg_04(end),Q_04(end,:),'DisplayName','a = 0.4');
plot(params.sigma_elem.*xg_06(end),Q_06(end,:),'DisplayName','a = 0.6');
plot(params.sigma_elem.*xg_08(end),Q_08(end,:),'DisplayName','a = 0.8');
plot(params.sigma_elem.*xg_1(end),Q_1(end,:),'DisplayName','a = 1');
plot(params.sigma_elem.*xg_12(end),Q_12(end,:),'DisplayName','a = 1.2');
plot(params.sigma_elem.*xg_14(end),Q_14(end,:),'DisplayName','a = 1.4');
plot(params.sigma_elem.*xg_16(end),Q_16(end,:),'DisplayName','a = 1.6');
plot(params.sigma_elem.*xg_18(end),Q_18(end,:),'DisplayName','a = 1.8');
plot(params.sigma_elem.*xg_2(end),Q_2(end,:),'DisplayName','a = 2');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("flow (m^3/s)");

subplot(3,1,2);
title("Effective pressures");
plot(params.sigma_elem.*xg_0(end),N_0(end,:),'DisplayName','a = 0');hold on;
plot(params.sigma_elem.*xg_02(end),N_02(end,:),'DisplayName','a = 0.2');
plot(params.sigma_elem.*xg_04(end),N_04(end,:),'DisplayName','a = 0.4');
plot(params.sigma_elem.*xg_06(end),N_06(end,:),'DisplayName','a = 0.6');
plot(params.sigma_elem.*xg_08(end),N_08(end,:),'DisplayName','a = 0.8');
plot(params.sigma_elem.*xg_1(end),N_1(end,:),'DisplayName','a = 1');
plot(params.sigma_elem.*xg_12(end),N_12(end,:),'DisplayName','a = 1.2');
plot(params.sigma_elem.*xg_14(end),N_14(end,:),'DisplayName','a = 1.4');
plot(params.sigma_elem.*xg_16(end),N_16(end,:),'DisplayName','a = 1.6');
plot(params.sigma_elem.*xg_18(end),N_18(end,:),'DisplayName','a = 1.8');
plot(params.sigma_elem.*xg_2(end),N_2(end,:),'DisplayName','a = 2');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("effective pressure (Pa)");


subplot(3,1,3);
title("Surface areas");
plot(params.sigma_elem.*xg_0(end),S_0(end,:),'DisplayName','a = 0');hold on;
plot(params.sigma_elem.*xg_02(end),S_02(end,:),'DisplayName','a = 0.2');
plot(params.sigma_elem.*xg_04(end),S_04(end,:),'DisplayName','a = 0.4');
plot(params.sigma_elem.*xg_06(end),S_06(end,:),'DisplayName','a = 0.6');
plot(params.sigma_elem.*xg_08(end),S_08(end,:),'DisplayName','a = 0.8');
plot(params.sigma_elem.*xg_1(end),S_1(end,:),'DisplayName','a = 1');
plot(params.sigma_elem.*xg_12(end),S_12(end,:),'DisplayName','a = 1.2');
plot(params.sigma_elem.*xg_14(end),S_14(end,:),'DisplayName','a = 1.4');
plot(params.sigma_elem.*xg_16(end),S_16(end,:),'DisplayName','a = 1.6');
plot(params.sigma_elem.*xg_18(end),S_18(end,:),'DisplayName','a = 1.8');
plot(params.sigma_elem.*xg_2(end),S_2(end,:),'DisplayName','a = 2');
legend('Location','eastoutside');
xlabel('distance (km)');ylabel("area (m^2)");


