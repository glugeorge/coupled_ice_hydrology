clear; close all;

load S1C_sensitivity_highres.mat
i = 4;
j = 4;
k = 3;
l = 4;
Q = ones(params.Nh,1);
N = ones(params.Nh,1);
S = ones(params.Nh,1); 
h = 2.*ones(length(params.sigma_elem),1);
xg = 290e3/params.x0; % Set high past sill for retreat
u = 0.1*(params.sigma_elem.^(1/3))+0.01; % 0.1 for C = 0.5, 0.3 for C = 0.1-0.4


init = [Q;N;S;h;u;xg];
[N_pos,max_h,xg,next_init,exitflag] = solve_steady_state(a_vals(i),C_vals(j),A_vals(k),M_vals(l),params,init);
N_pos_arr(i,j,k,l) = N_pos;
max_h_arr(i,j,k,l) = max_h;
xg_arr(i,j,k,l) = xg;

%% Plotting
fig = figure(1);
t = tiledlayout(length(a_vals),length(C_vals),'TileSpacing','Compact');
for i=1:length(a_vals)
    for j=1:length(C_vals)
        nexttile;
        heatmap(A_vals,M_vals,squeeze(N_pos_arr(i,j,:,:))','ColorbarVisible','off');
        %set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %clim([0.9,1]);
        title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])
    end

end
%cb = colorbar(); 
%cb.Layout.Tile = 'east'; % Assign colorbar location
%%
fig = figure(2);
t = tiledlayout(length(a_vals),length(C_vals),'TileSpacing','Compact');
for i=1:length(a_vals)
    for j=1:length(C_vals)
        nexttile;
        heatmap(A_vals,M_vals,squeeze(xg_arr(i,j,:,:))','ColorbarVisible','off');
        
        %clim([0.95,1]);
        title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])
    end

end
%% Save values
save("S1C_sensitivity_highres.mat");


%% Solve function
function [N_pos,max_h,xg,next_init,exitflag] = solve_steady_state(a,C,A,M,params_init,init)
bed_func = @bed_flat;
params = params_init;
params.b0 = -100;           %bed topo at x=0
params.bx = -1e-3;          %linear bed slope

params.sill_min = 2000e3;   %sill min x position
params.sill_max = 2100e3;   %sill max x position
params.sill_slope = 1e-3;   %slope of sill
%% Physical parameters
params.Cf = 1;
params.A = A;
params.n = 3;
params.rho_i = 917;
params.rho_w = 1028;
params.g = 9.81;
params.C = C; 
params.As = 2.26e-21;
params.f = 0.07; % From Kingslake thesis
params.K0 = 10^-24; % From Kingslake thesis  
params.L = 3.3e5; % Kingslake thesis
params.year = 3600*24*365;
%% Scaling params (coupled model equations solved in non-dim form)

params.psi0 = params.rho_w*params.g*params.h0/params.x0;
params.M0 = params.Q0/params.x0;
params.m0 = params.Q0*params.psi0/params.L;
params.eps_r = params.m0*params.x0/(params.rho_i*params.Q0);
params.S0 = (params.f*params.rho_w*params.g*params.Q0^2/params.psi0)^(3/8);
params.th0 = params.rho_i*params.S0/params.m0;
params.N0 = (params.K0*params.th0)^(-1/3);
params.delta = params.N0/(params.x0*params.psi0);
%params.u0 = (rho_i*g*params.h0^2/(C*params.N0*params.x0))^n;
params.u0 = params.As*(params.C*params.N0)^params.n;
params.t0 = params.x0/params.u0;
params.a0 = params.h0/params.t0;
params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));
%gamma = As*(C*N0)^n;
params.gamma = (params.C*params.N0*params.x0)/(params.rho_i*params.g*params.h0^2);
params.beta = params.th0/params.t0;
params.r = params.rho_i/params.rho_w;

params.transient = 0;

%% Establish timings
params.year = 3600*24*365;  %number of seconds in a year
params.Nt =10;                    %number of time steps - normally 150
params.end_year = 10; %normally 7500

params.dt = params.end_year*params.year/params.Nt;

%% Determine at what points there is coupling
% 1 - coupling on, 0 - coupling off
params.hydro_u_from_ice_u = 1;
params.hydro_psi_from_ice_h = 1;
params.ice_N_from_hydro = 1;

%% Initial "steady state" conditions
params.shear_scale = 1;

params.M = M/params.M0; % zero when using schoof bed
params.N_terminus = 0;
params.accum = a./params.year;
%hf = (-bed_func(xg.*params.x0,params)/params.h0)/params.r;
params.Q_in = 0.001/params.Q0;


params.ice_start = 3*params.Nh;


QNShuxg0 = init;
Q = QNShuxg0(1:params.Nh);
N = QNShuxg0(params.Nh+1:2*params.Nh);
S = QNShuxg0(2*params.Nh+1:3*params.Nh);
h = QNShuxg0(params.ice_start+1:params.ice_start+ params.Nx);
u = QNShuxg0(params.ice_start + params.Nx+1:params.ice_start+2*params.Nx);
xg = QNShuxg0(params.ice_start+2*params.Nx+1);
params.S_old = S;
params.h_old = h;
params.xg_old = xg;
sig_old = params.sigma;
sige_old = params.sigma_elem;

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
flf = @(QNShuxg) ice_hydro_equations_coulomb(QNShuxg,params,bed_func);

[QNShuxg_init,F,exitflag,output,JAC] = fsolve(flf,QNShuxg0,options);

Q = QNShuxg_init(1:params.Nh);
N = QNShuxg_init(params.Nh+1:2*params.Nh);
S = QNShuxg_init(2*params.Nh+1:3*params.Nh);
h = QNShuxg_init(params.ice_start+1:params.ice_start+ params.Nx);
u = QNShuxg_init(params.ice_start + params.Nx+1:params.ice_start+2*params.Nx);
xg = QNShuxg_init(params.ice_start+2*params.Nx+1);
hf = (-bed_func(xg.*params.x0,params)/params.h0)/(params.r);

[N_max,I] = max(N);
N_pos = params.sigma_h(I);
max_h = max(h);
next_init = QNShuxg_init;

end