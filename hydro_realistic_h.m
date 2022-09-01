% Solving hydrology equations using finite differences
clear all;
% Using new scales + realistic h
% I think this would be the initial model to look at as a baseline

% Load in determined steady state h, u 
load init_cond.mat

h_grid = params.sigma_elem*xg; % nondimensionalized xgrid for h
u_grid = params.sigma*xg; % nondimensionalized xgrid for u
params.xg = xg;
params.Nt = 1000;                    %number of time steps
params.year = 3600*24*365;  %number of seconds in a year
params.dt = 0.05;
params.Nx = 200;                    %number of grid points
%params.N1 = 100;                    %number of grid points in coarse domain
params.sigGZ = 1; % 0.97;                %extent of coarse grid (where GL is at sigma=1)
%sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
%sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
%params.sigma = [sigma1, sigma2(2:end)]';  
params.sigma = linspace(0,1,params.Nx)';%[linspace(params.sigGZ/(params.Nx+0.5),params.sigGZ,params.Nx)]';
params.dsigma = diff(params.sigma); %grid spacing

% interpolate h to same grid
h_interp = interp1(h_grid,h,params.sigma.*xg,'linear','extrap');
% scale constants
% Define given constants
A = 4.227e-25; % From Alex Robel's code
n = 3;
rho_i = 917;
rho_w = 1028;
g = 9.81;
C = 0.84*0.5; % From Hewitt and Fowler, 2008
As = 2*2.4e-24/(0.5*C^n); % Calculated 
p = 1/3; % From Hewitt's Karthaus slides 
q = 1/3; % From Hewitt's Karthaus slides
f = 0.07; % From Kingslake thesis
K0 = 10^-24; % From Kingslake thesis  
L = 3.3e5; % Kingslake thesis


params.x0 = 10*10^3;
params.h0 = 100;
params.Q0 = 1500;

params.psi0 = rho_w*g*0.001;
params.M0 = params.Q0/params.x0;
params.m0 = params.Q0*params.psi0/L;
params.eps_r = params.m0*params.x0/(rho_i*params.Q0);
params.S0 = (f*rho_w*g*params.Q0^2/params.psi0)^(3/8);
params.th0 = rho_i*params.S0/params.m0;
params.N0 = (K0*params.th0)^(-1/3);
params.delta = params.N0/(params.x0*params.psi0);
%params.u0 = (rho_i*g*params.h0^2/(C*params.N0*params.x0))^n;
params.u0 = As*(C*params.N0)^n;
params.t0 = params.x0/params.u0;
params.a0 = params.h0/params.t0;
params.alpha = 2*params.u0^(1/n)/(rho_i*g*params.h0*(params.x0*A)^(1/n));
%gamma = As*(C*N0)^n;
params.gamma = (C*params.N0*params.x0)/(rho_i*g*params.h0^2);
params.beta = params.th0/params.t0;

% Initial conditions
% figure out initial N
Q = 0.001*ones(params.Nx,1);
N = ones(params.Nx,1);
S = 5/params.S0*ones(params.Nx,1); 
params.S_old = S;
params.u = params.th0*0.*ones(params.Nx,1)./(params.year*params.x0); % Set 0 velocity first
% define psi with realistic h
phi_b = 0.001; % slope
p = rho_i*g*h_interp*params.h0;
params.psi = (rho_w*g*sin(phi_b)-gradient(p)./gradient(params.sigma.*params.x0.*xg))./params.psi0;

params.M = 5*10^-4/params.M0; 
params.N_terminus = 0;%rho_i*g*h_interp(end)*params.h0/params.N0;
params.r = rho_i/rho_w;

QNS = [Q; N; S];

Qs = nan.*ones(params.Nt,params.Nx);
Ns = nan.*ones(params.Nt,params.Nx);
Ss = nan.*ones(params.Nt,params.Nx);

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);

for t=2:params.Nt
    t
    eqns = @(QNS) hydro_eqns(QNS,params);
    [QNS,F,exitflag,output,JAC] = fsolve(eqns,QNS,options);
    Q = QNS(1:params.Nx);
    N = QNS(params.Nx+1:2*params.Nx);
    S = QNS(2*params.Nx+1:3*params.Nx);
    params.S_old = S;
    Qs(t,:)=Q;
    Ns(t,:)=N;
    Ss(t,:)=S;
end


%% Plotting stuff
Q_dim = Qs.*params.Q0;
N_dim = Ns.*params.N0;
S_dim = Ss.*params.S0;
ts = linspace(0,params.Nt*params.dt*params.th0./params.year,params.Nt);
disp(xg*params.x0);
figure()
subplot(3,1,1);surface(ts,params.sigma,Q_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Flow (m^3/s)');set(gca,'Ydir','Reverse')
subplot(3,1,2);surface(ts,params.sigma,N_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Effective Pressure (Pa)');set(gca,'Ydir','Reverse')
subplot(3,1,3);surface(ts,params.sigma,S_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Surface Area (m^2)');set(gca,'Ydir','Reverse')

figure('Name','Evolution of Drainage');
subplot(3,1,1);
plot(ts,N_dim(:,1),'DisplayName','Channel Entrance');
hold on;
plot(ts,N_dim(:,end),'--','DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)','Interpreter','latex');
ylabel('Effective Pressure, \emph{N} (Pa)','Interpreter','latex');
title('Effective Pressure over time');

subplot(3,1,2);
plot(ts,Q_dim(:,1),'DisplayName','Channel Entrance');
hold on;
plot(ts,Q_dim(:,end),'DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)','Interpreter','latex');
ylabel('channel flow rate, \emph{Q} $(m^{3} s^{-1})$','Interpreter','latex');
title('Flow over time');
legend('Location','northwest');

subplot(3,1,3);
plot(ts,S_dim(:,1),'-o','DisplayName','Channel Entrance');
hold on;
plot(ts,S_dim(:,end),'-o','DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)','Interpreter','latex');
ylabel('channel cross-sectional area, \emph{S} $(m^2)$','Interpreter','latex');
title('Channel area over time');
legend('Location','northwest');

save 'realistic_h.mat' Q_dim N_dim S_dim ts

%% Function to solve
function F = hydro_eqns(QNS,params)
    M = params.M;
    Q_in = 10/params.Q0;
    Nx = params.Nx;
    Q = QNS(1:Nx);
    N = QNS(Nx+1:2*Nx);
    S = QNS(2*Nx+1:3*Nx);
    xg = params.xg;
    xg_old = xg;
    sigma = params.sigma;
    % Q
    fq(1) = Q(1) - Q_in; % Boundary condition

    fq(2:Nx-1) = (params.eps_r*(params.r-1))*abs(Q(2:Nx-1)).^3./S(2:Nx-1).^(8/3) +...
                        params.eps_r.*S(2:Nx-1).*N(2:Nx-1).^3 + M + ...
                        params.eps_r.*params.beta.*params.u(2:Nx-1).*(S(3:Nx)-S(1:Nx-2))./(2*xg*params.dsigma(2:Nx-1)) - ...
                        (Q(3:Nx)-Q(1:Nx-2))./(2*xg*params.dsigma(2:Nx-1));
    fq(Nx) = (params.eps_r*(params.r-1))*abs(Q(Nx)).^3./S(Nx).^(8/3) +...
                        params.eps_r*S(Nx).*N(Nx).^3 + M + ...
                        params.eps_r.*params.beta.*params.u(Nx).*(S(Nx)-S(Nx-1))./(xg*params.dsigma(Nx-1)) - ...
                        (Q(Nx)-Q(Nx-1))./(xg*params.dsigma(Nx-1)); % one sided difference instead
    % N 
    fn(1) = Q(1).*abs(Q(1))./(S(1).^(8/3)) - params.psi(1) - ...
                        params.delta*(N(2)-N(1))./(xg*params.dsigma(1)); % use 1 sided difference instead of symmetry argument
    fn(2:Nx-1) =  Q(2:Nx-1).*abs(Q(2:Nx-1))./(S(2:Nx-1).^(8/3)) - params.psi(2:Nx-1) - ...
                        params.delta*(N(3:Nx)-N(1:Nx-2))./(2*xg*params.dsigma(2:Nx-1));
    fn(Nx) = N(Nx) - params.N_terminus; % Boundary condition

    % S
    fs(1) = abs(Q(1)).^3./(S(1).^(8/3)) - ... 
                        S(1).*N(1).^3 + ...
                        (sigma(1)*(xg-xg_old)/params.dt - params.beta.*params.u(1)).*(S(2)-S(1))./(xg*params.dsigma(1)) - ...
                        (S(1)-params.S_old(1))./params.dt; % one sided differece
    fs(2:Nx-1)= abs(Q(2:Nx-1)).^3./(S(2:Nx-1).^(8/3)) - ... 
                        S(2:Nx-1).*N(2:Nx-1).^3 + ...
                        (sigma(2:Nx-1).*(xg-xg_old)./params.dt - params.beta.*params.u(2:Nx-1)).*(S(3:Nx)-S(1:Nx-2))./(2*xg*params.dsigma(2:Nx-1)) - ...
                        (S(2:Nx-1)-params.S_old(2:Nx-1))./params.dt; 
    fs(Nx)= abs(Q(Nx)).^3./(S(Nx).^(8/3)) - ... 
                        S(Nx).*N(Nx).^3 + ...
                        (sigma(Nx).*(xg-xg_old)./params.dt - params.beta.*params.u(Nx)).*(S(Nx)-S(Nx-1))./(xg*params.dsigma(Nx-1)) - ...
                        (S(Nx)-params.S_old(Nx))./params.dt; % one sided differece

    F = [fq;fn;fs];
end
