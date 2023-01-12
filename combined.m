close all
clear all
clc

%% Bed parameters
params.b0 = -100;           %bed topo at x=0
params.bx = -1e-3;          %linear bed slope

params.sill_min = 2000e3;   %sill min x position
params.sill_max = 2100e3;   %sill max x position
params.sill_slope = 1e-3;   %slope of sill
  
%% Physical parameters
params.A = 2.9e-25; 
params.n = 3;
params.rho_i = 917;
params.rho_w = 1028;
params.g = 9.81;
params.C = 0.2; % base 0.2
params.As = 2.26e-21; % Calculated 
params.f = 0.07; % From Kingslake thesis
params.K0 = 10^-24; % From Kingslake thesis  
params.L = 3.3e5; % Kingslake thesis
params.year = 3600*24*365;
%% Scaling params (coupled model equations solved in non-dim form)
params.x0 = 100*10^3;
params.h0 = 100;
params.Q0 = 1500;

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

%% Grid parameters - ice sheet
params.Nx = 1002;                    %number of grid points - 200
params.N1 = 2;                    %number of grid points in coarse domain - 100
params.Nh = 1002;
params.sigGZ = 0.01;                %extent of coarse grid (where GL is at sigma=1) - 0.97
sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
params.sigma = [sigma1, sigma2(2:end)]';    %grid points on velocity (includes GL, not ice divide)
params.dsigma = diff(params.sigma); %grid spacing
params.sigma_elem = [0;(params.sigma(1:params.Nx-1) + params.sigma(2:params.Nx))./2 ]; %grid points on thickness (includes divide, not GL)

%% Grid parameters - hydro
params.sigma_h = linspace(0,1,params.Nh)';
params.dsigma_h = diff(params.sigma_h); %grid spacing

%% Establish timings
params.year = 3600*24*365;  %number of seconds in a year
params.Nt =50;                    %number of time steps - normally 150
params.end_year = 1000; %normally 7500

params.dt = params.end_year*params.year/params.Nt;

%% Determine at what points there is coupling
% 1 - coupling on, 0 - coupling off
params.hydro_u_from_ice_u = 1;
params.hydro_psi_from_ice_h = 1;
params.ice_N_from_hydro = 1;

%% Initial "steady state" conditions
params.shear_scale = 0.1;
Q = 0.001*ones(params.Nh,1);
N = ones(params.Nh,1);
S = 5/params.S0*ones(params.Nh,1); 
params.S_old = S;
params.M = 0e-4/params.M0; % zero when using schoof bed
params.N_terminus = 0;
params.accum = 1./params.year;
xg = 100e3/params.x0;
hf = (-bed(xg.*params.x0,params)/params.h0)/params.r;
h =  1 - (1-hf).*params.sigma;
u = 0.3*(params.sigma_elem.^(1/3)) + 1e-3; % 0.1 for C = 0.5, 0.3 for C = 0.1-0.4
params.Q_in = 10/params.Q0;

params.h_old = h;
params.xg_old = xg;
params.ice_start = 3*params.Nh;

sig_old = params.sigma;
sige_old = params.sigma_elem;
QNShuxg0 = [Q;N;S;h;u;xg];

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
flf = @(QNShuxg) combined_hydro_ice_eqns(QNShuxg,params);

[QNShuxg_init,F,exitflag,output,JAC] = fsolve(flf,QNShuxg0,options);

Q = QNShuxg_init(1:params.Nh);
N = QNShuxg_init(params.Nh+1:2*params.Nh);
S = QNShuxg_init(2*params.Nh+1:3*params.Nh);
h = QNShuxg_init(params.ice_start+1:params.ice_start+ params.Nx);
u = QNShuxg_init(params.ice_start + params.Nx+1:params.ice_start+2*params.Nx);
xg = QNShuxg_init(params.ice_start+2*params.Nx+1);
hf = (-bed(xg.*params.x0,params)/params.h0)/(params.r);

%% Final steady state solution
% params.accum 
% = 1./params.year;
% params.Q_in = 10/params.Q0;
 params.A = 4.9e-25; 
 params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));
 flf = @(QNShuxg) combined_hydro_ice_eqns(QNShuxg,params);
 [QNShuxg_final,F,exitflag,output,JAC] = fsolve(flf,QNShuxg_init,options);
 xg_f = QNShuxg_final(params.ice_start+2*params.Nx+1);

%% Now for evolution 
Qs = nan.*ones(params.Nt,params.Nh);
Ns = nan.*ones(params.Nt,params.Nh);
Ss = nan.*ones(params.Nt,params.Nh);
hs = nan.*ones(params.Nt,params.Nx);
us = nan.*ones(params.Nt,params.Nx);
xgs = nan.*ones(1,params.Nt);
QNShuxg_t = QNShuxg_init;

Qs(1,:) = QNShuxg_t(1:params.Nh);
Ns(1,:) = QNShuxg_t(params.Nh+1:2*params.Nh);
Ss(1,:) = QNShuxg_t(2*params.Nh+1:3*params.Nh);
hs(1,:) = QNShuxg_t(params.ice_start+1:params.ice_start+ params.Nx);
us(1,:) = QNShuxg_t(params.ice_start+params.Nx+1:params.ice_start+2*params.Nx);
xgs(1) = QNShuxg_t(params.ice_start+2*params.Nx+1);

params.h_old = h;
params.xg_old =xg;
params.S_old = S;
params.transient = 1;
time_to_ss = 0; % to 99 percent
for t=2:params.Nt
%     if t == 50
%         params.A = 0.9e-25; 
%         params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));
% 
%     end
    flf = @(QNShuxg) combined_hydro_ice_eqns(QNShuxg,params);
    [QNShuxg_t,F,exitflag,output,JAC] = fsolve(flf,QNShuxg_t,options);
    t
    Qs(t,:) = QNShuxg_t(1:params.Nh);
    Ns(t,:) = QNShuxg_t(params.Nh+1:2*params.Nh);
    Ss(t,:) = QNShuxg_t(2*params.Nh+1:3*params.Nh);
    hs(t,:) = QNShuxg_t(params.ice_start+1:params.ice_start+ params.Nx);
    us(t,:) = QNShuxg_t(params.ice_start+params.Nx+1:params.ice_start+2*params.Nx);
    xgs(t) = QNShuxg_t(params.ice_start+2*params.Nx+1);
    params.h_old = QNShuxg_t(params.ice_start+1:params.ice_start+ params.Nx);
    params.xg_old = xgs(t);
    params.S_old = QNShuxg_t(2*params.Nh+1:3*params.Nh);
    %if abs(xg_f - xgs(t)) < 0.001*xg_f && time_to_ss == 0
    %    time_to_ss = (t-1)*params.dt/params.year;
    %end
end

%% Retreat (back over)
% params.A = 2.9e-25; 
% params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));

%% Plotting
ts = linspace(0,params.end_year,params.Nt);
figure();
subplot(3,1,1);plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time (yr)');ylabel('x_g');
subplot(3,1,2);contourf(ts,params.sigma_elem,hs'.*params.h0);colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse');
subplot(3,1,3);contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse');

figure();
Q_dim = Qs.*params.Q0;
N_dim = Ns.*params.N0;
S_dim = Ss.*params.S0;
subplot(3,1,1);surface(ts,params.sigma_h,Q_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Flow (m^3/s)');set(gca,'Ydir','Reverse')
subplot(3,1,2);surface(ts,params.sigma_h,N_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Effective Pressure (Pa)');set(gca,'Ydir','Reverse')
subplot(3,1,3);surface(ts,params.sigma_h,S_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Surface Area (m^2)');set(gca,'Ydir','Reverse')

%% Saving values
results.params = params;
results.init_cond = QNShuxg_init;
%results.steady_state = QNShuxg_final;
results.xgs = xgs;
results.ts = ts;
results.hs = hs';
results.us = us';
results.Qs = Qs';
results.Ns = Ns';
results.Ss = Ss';
%results.time_to_ss = time_to_ss; 

%fname = strcat('base_run',num2str(params.A*1e25),'_c.mat');
fname = strcat('C_',num2str(params.C),'_A_',num2str(params.A*1e25),'_c_0.1shear.mat');
%fname = strcat('Nh_',num2str(params.Nh),'_coarse_',num2str(params.N1),'_fine_',num2str(params.Nx-params.N1),'.mat');
save(fname,'results');