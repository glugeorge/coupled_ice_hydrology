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
params.Q0 = 10;

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
params.Nx = 1000;                    %number of grid points - 200
params.N1 = 100;                    %number of grid points in coarse domain - 100
params.Nh = 1000;
params.sigGZ = 0.73;                %extent of coarse grid (where GL is at sigma=1) - 0.97
sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
params.sigma = [sigma1, sigma2(2:end)]';    %grid points on velocity (includes GL, not ice divide)
params.dsigma = diff(params.sigma); %grid spacing
params.sigma_elem = [0;(params.sigma(1:params.Nx-1) + params.sigma(2:params.Nx))./2 ]; %grid points on thickness (includes divide, not GL)

%% Grid parameters - hydro
params.sigma_h = linspace(0,1,params.Nh)';
params.dsigma_h = diff(params.sigma_h); %grid spacing


%% Define velocity and ice shape
% load C_0.2_A_4.9_c_nosplit.mat
% params = results.params;
% params.xg = results.xgs(end);
% params.h = results.hs(:,end);
% params.u = results.us(:,end);
% %params.xg = 100e3./params.x0;
params.xg = 100e3./params.x0;

params.h =  10*sqrt(1-params.sigma_elem); %ones(params.Nx,1).*1000./params.h0;
%params.h(params.N1+1) = (1+0.5*(1+(params.dsigma(params.N1)/params.dsigma(params.N1-1))))*params.h(params.N1) - 0.5*(1+(params.dsigma(params.N1)/params.dsigma(params.N1-1)))*params.h(params.N1-1);
params.u = (500/params.year)*ones(params.Nx,1)./params.u0;

%% Establish timings
params.year = 3600*24*365;  %number of seconds in a year
params.Nt =1;                    %number of time steps - normally 150
params.end_year = 30; %normally 7500

params.dt = params.end_year*params.year/params.Nt;

%% Initial condition
Q = 10*ones(params.Nh,1)./params.Q0;
N = ones(params.Nh,1);
S = 5/params.S0*ones(params.Nh,1); 
params.S_old = S;
params.M = 0e-4/params.M0; % zero when using schoof bed
params.N_terminus = 0;
params.accum = 1./params.year;

params.Q_in = 10/params.Q0;

params.ice_start = 3*params.Nh;

QNShuxg_init = [Q;N;S];%;h];u;xg];

Qs = nan.*ones(params.Nt,params.Nh);
Ns = nan.*ones(params.Nt,params.Nh);
Ss = nan.*ones(params.Nt,params.Nh);
% hs = nan.*ones(params.Nt,params.Nx);
% us = nan.*ones(params.Nt,params.Nx);
%xgs = nan.*ones(1,params.Nt);
QNShuxg_t = QNShuxg_init;

Qs(1,:) = QNShuxg_t(1:params.Nh);
Ns(1,:) = QNShuxg_t(params.Nh+1:2*params.Nh);
Ss(1,:) = QNShuxg_t(2*params.Nh+1:3*params.Nh);
%hs(1,:) = QNShuxg_t(params.ice_start+1:params.ice_start+ params.Nx);
%us(1,:) = QNShuxg_t(params.ice_start+params.Nx+1:params.ice_start+2*params.Nx);
%xgs(1) = QNShuxg_t(params.ice_start+2*params.Nx+1);

params.h_old = params.h;
params.xg_old =params.xg;
params.S_old = S;
%params.transient = 1;
options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);

% Evolve
for t=1:params.Nt

    flf = @(QNShuxg) hydro_eqns(QNShuxg,params);
    [QNShuxg_t,F,exitflag,output,JAC] = fsolve(flf,QNShuxg_t,options);
    t
    Qs(t,:) = QNShuxg_t(1:params.Nh);
    Ns(t,:) = QNShuxg_t(params.Nh+1:2*params.Nh);
    Ss(t,:) = QNShuxg_t(2*params.Nh+1:3*params.Nh);
    %hs(t,:) = QNShuxg_t(params.ice_start+1:params.ice_start+ params.Nx);
    %us(t,:) = QNShuxg_t(params.ice_start+params.Nx+1:params.ice_start+2*params.Nx);
    %xgs(t) = QNShuxg_t(params.ice_start+2*params.Nx+1);
    params.h_old = params.h;%QNShuxg_t(params.ice_start+1:params.ice_start+ params.Nx);
    params.xg_old = params.xg; %xgs(t);
    params.S_old = QNShuxg_t(2*params.Nh+1:3*params.Nh);

end

%% Plotting
ts = linspace(0,params.end_year,params.Nt);

figure();
Q_dim = Qs.*params.Q0;
N_dim = Ns.*params.N0;
S_dim = Ss.*params.S0;
subplot(3,1,1);surface(ts,params.sigma_h,Q_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Flow (m^3/s)');set(gca,'Ydir','Reverse')
subplot(3,1,2);surface(ts,params.sigma_h,N_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Effective Pressure (Pa)');set(gca,'Ydir','Reverse')
subplot(3,1,3);surface(ts,params.sigma_h,S_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Surface Area (m^2)');set(gca,'Ydir','Reverse')

figure();
ax1 = subplot(5,1,1); 
plot(params.sigma_h,Qs(end,:)); ylabel('Q');title('Nondimensionalized steady state variables');
hold on;
ax2 =subplot(5,1,2);
plot(params.sigma_h,Ns(end,:));ylabel('N');hold on;
ax3 = subplot(5,1,3);
plot(params.sigma_h,Ss(end,:));ylabel('S');hold on;
ax4 = subplot(5,1,4);
plot(params.sigma_elem,params.h);ylabel('h');hold on;
plot(linspace(0,1),bed(linspace(0,1).*params.xg,params)./params.h0,'-k'); 
ax5 = subplot(5,1,5);
plot(params.sigma,params.u);hold on;ylabel('u');xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
linkaxes([ax1,ax2,ax3,ax4,ax5],'x')

%% Saving relevant variables
save quad_oneway_1.mat params Qs Ns Ss;

%% Functions
function F = hydro_eqns(QNShuxg,params)
    % unpack variables
    M = params.M;
    Q_in = params.Q_in;
    Nx = params.Nx;
    Nh = params.Nh;
    ice_start = params.ice_start;
    Q = QNShuxg(1:Nh);
    N = QNShuxg(Nh+1:2*Nh);
    S = QNShuxg(2*Nh+1:3*Nh);
    h = params.h; %QNShuxg(ice_start+1:ice_start + Nx);
    u = params.u; %QNShuxg(ice_start + Nx+1:ice_start + 2*Nx);
    xg = params.xg; %QNShuxg(ice_start + 2*Nx+1);
    hf = (-bed(xg.*params.x0,params)/params.h0)/(params.r);
    xg_old = xg;
    %grid params unpack
    dt = params.dt/params.t0;
    dth= params.dt/params.th0;
    ds = params.dsigma;
    Nx = params.Nx;
    N1 = params.N1;
    sigma = params.sigma;
    sigma_elem = params.sigma_elem;
    b = -bed(xg.*sigma.*params.x0,params)/params.h0;
    dbdx = gradient(bed(xg.*params.sigma_h.*params.x0,params))./gradient(params.sigma_h.*params.x0*xg);
    Fh = zeros(Nx,1);
    Fu = zeros(Nx,1);
    Fxg = 0; 
    %physical params unpack
    m     = 1/params.n;
    nglen = params.n;
    accum = params.accum;
    a = accum/params.a0;
    gamma = params.gamma;
    alpha = params.alpha;
    ss = params.transient;
       
    % Assign values based off coupling parameters
    u_ice_interp = interp1(sigma,u,params.sigma_h,"linear","extrap");

    h_interp = interp1(sigma_elem,h.*params.h0,params.sigma_h,'linear','extrap');
    params.psi = (-params.rho_w*params.g.*dbdx-params.rho_i*params.g.*gradient(h_interp)./gradient(params.sigma_h.*params.x0*xg))./params.psi0;
    

    % Q
    fq(1) = Q(1) - Q_in; % Boundary condition

    fq(2:Nh-1) = (params.eps_r*(params.r-1))*abs(Q(2:Nh-1)).^3./S(2:Nh-1).^(8/3) +...
                        params.eps_r.*S(2:Nh-1).*N(2:Nh-1).^3 + M + ...
                        params.eps_r.*params.beta.*u_ice_interp(2:Nh-1).*(S(3:Nh)-S(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1)) - ...
                        (Q(3:Nh)-Q(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1));
    fq(Nh) = (params.eps_r*(params.r-1))*abs(Q(Nh)).^3./S(Nh).^(8/3) +...
                        params.eps_r*S(Nh).*N(Nh).^3 + M + ...
                        params.eps_r.*params.beta.*u_ice_interp(Nh).*(S(Nh)-S(Nh-1))./(xg*params.dsigma_h(Nh-1)) - ...
                        (Q(Nh)-Q(Nh-1))./(xg*params.dsigma_h(Nh-1)); % one sided difference instead
    % N 
    fn(1) = Q(1).*abs(Q(1))./(S(1).^(8/3)) - params.psi(1) - ...
                        params.delta*(N(2)-N(1))./(xg*params.dsigma_h(1)); % use 1 sided difference instead of symmetry argument
    fn(2:Nh-1) =  Q(2:Nh-1).*abs(Q(2:Nh-1))./(S(2:Nh-1).^(8/3)) - params.psi(2:Nh-1) - ...
                        params.delta*(N(3:Nh)-N(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1));
    fn(Nh) = N(Nh) - params.N_terminus; % Boundary condition
    
    % S
    fs(1) = abs(Q(1)).^3./(S(1).^(8/3)) - ... 
                        S(1).*N(1).^3 + ...
                        (ss.*params.sigma_h(1)*(xg-xg_old)/dth - params.beta.*u_ice_interp(1)).*(S(2)-S(1))./(xg*params.dsigma_h(1)) - ...
                        ss.*(S(1)-params.S_old(1))./dth; % one sided difference
    fs(2:Nh-1)= abs(Q(2:Nh-1)).^3./(S(2:Nh-1).^(8/3)) - ... 
                        S(2:Nh-1).*N(2:Nh-1).^3 + ...
                        (ss.*params.sigma_h(2:Nh-1).*(xg-xg_old)./dth - params.beta.*u_ice_interp(2:Nh-1)).*(S(3:Nh)-S(1:Nh-2))./(2*xg*params.dsigma_h(2:Nh-1)) - ...
                        ss.*(S(2:Nh-1)- params.S_old(2:Nh-1))./dth; 
    fs(Nh)= abs(Q(Nh)).^3./(S(Nh).^(8/3)) - ... 
                        S(Nh).*N(Nh).^3 + ...
                        (ss.*params.sigma_h(Nh).*(xg-xg_old)./dth - params.beta.*u_ice_interp(Nh)).*(S(Nh)-S(Nh-1))./(xg*params.dsigma_h(Nh-1)) - ...
                        ss.*(S(Nh)-params.S_old(Nh))./dth; % one sided difference

    F = [fq';fn';fs'];
end