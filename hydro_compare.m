%% Establish consistent time scales
endTime = 10*60*60*24*365.25; % five years
params.Nt = 5000;                    %number of time steps
params_new.Nt = 5000;                    %number of time steps
params.Nx = 200;  
params_new.Nx = 200;

%% Old Using scales from Kingslake et al to compare with BVP solver

%params.N1 = 100;                    %number of grid points in coarse domain
params.year = 3600*24*365;  %number of seconds in a year
params.sigGZ = 1; % 0.97;                %extent of coarse grid (where GL is at sigma=1)
%sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
%sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
%params.sigma = [sigma1, sigma2(2:end)]';  
params.sigma = linspace(0,1,params.Nx)';%[linspace(params.sigGZ/(params.Nx+0.5),params.sigGZ,params.Nx)]';
params.dsigma = diff(params.sigma); %grid spacing
% Thickness has this other sigma_elem used, need to dissect that

% scale constants
% Define given constants
rho_w = 1000; % kg/m^3 (water density)
rho_i = 900; %kg/m^3 (ice density)
g = 10; % m/s^2 (gravitational acceleration)
L = 3.3*10^5; % J/kg
K0 = 10^-24; % Pa^-3 s^-2
phi_s = 0.01; % surface slope
psi_0 = rho_w*g*phi_s;
n_prime = 0.1; % m^-1/3 s (hydraulic roughness)
f = 6.6*n_prime^2; 
s_to_y = 60*60*24*365.25;


params.x0 = 10*10^3;
params.h0 = 100;
params.Q0 = 1500;

params.m0 = params.Q0*psi_0/L;
params.M0 = params.Q0/params.x0; 
params.S0 = (f*rho_w*g*params.Q0^2/psi_0)^(3/8);
params.th0 = rho_i*params.S0*L/(psi_0*params.Q0); 
params.N0 = (K0*rho_i*params.S0*L/(psi_0*params.Q0))^(-1/3);
params.eps_r = params.x0*params.m0/(params.Q0*rho_i);
params.r = rho_i/rho_w;
params.delta = params.N0/(params.x0*psi_0);
params.psi = 1;
params.M = 5*10^-4/params.M0; 
params.N_terminus = rho_i*g*1/params.N0;
% Initial conditions
% figure out initial N
Q = 0.001*ones(params.Nx,1);
N = ones(params.Nx,1);
S = 5/params.S0*ones(params.Nx,1); 
params.S_old = S;
params.u = params.th0*1000.*ones(params.Nx,1)./(params.year*params.x0); % Set 0 velocity first
params.psi = params.psi*(1-3*exp(-20.*params.sigma));
QNS = [Q; N; S];

Qs = nan.*ones(params.Nt,params.Nx);
Ns = nan.*ones(params.Nt,params.Nx);
Ss = nan.*ones(params.Nt,params.Nx);

params.dt = endTime/(params.Nt*params.th0);

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);

%% NEW 
params_new.sigGZ = 1; % 0.97;                %extent of coarse grid (where GL is at sigma=1)
params_new.year = 3600*24*365;  %number of seconds in a year
%sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
%sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
%params.sigma = [sigma1, sigma2(2:end)]';  
params_new.sigma = linspace(0,1,params_new.Nx)';%[linspace(params.sigGZ/(params.Nx+0.5),params.sigGZ,params.Nx)]';
params_new.dsigma = diff(params_new.sigma); %grid spacing
% Thickness has this other sigma_elem used, need to dissect that

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


params_new.x0 = 10*10^3;
params_new.h0 = 100;
params_new.Q0 = 1500;

params_new.psi0 = rho_i*g*params_new.h0/params_new.x0;
params_new.M0 = params_new.Q0/params_new.x0;
params_new.m0 = params_new.Q0*params_new.psi0/L;
params_new.eps_r = params_new.m0*params_new.x0/(rho_i*params_new.Q0);
params_new.S0 = (f*rho_w*g*params_new.Q0^2/params_new.psi0)^(3/8);
params_new.th0 = rho_i*params_new.S0/params_new.m0;
params_new.N0 = (K0*params_new.th0)^(-1/3);
params_new.delta = params_new.N0/(params_new.x0*params_new.psi0);
%params_new.u0 = (rho_i*g*params_new.h0^2/(C*params_new.N0*params_new.x0))^n;
params_new.u0 = As*(C*params_new.N0)^n;
params_new.t0 = params_new.x0/params_new.u0;
params_new.a0 = params_new.h0/params_new.t0;
params_new.alpha = 2*params_new.u0^(1/n)/(rho_i*g*params_new.h0*(params_new.x0*A)^(1/n));
%gamma = As*(C*N0)^n;
params_new.gamma = (C*params_new.N0*params_new.x0)/(rho_i*g*params_new.h0^2);
params_new.beta = params_new.th0/params_new.t0;
% Initial conditions
% figure out initial N
Q = 0.001*ones(params_new.Nx,1);
N = ones(params_new.Nx,1);
S = 5/params_new.S0*ones(params_new.Nx,1); 
params_new.S_old = S;
params_new.u = params_new.th0*1000.*ones(params_new.Nx,1)./(params_new.year*params_new.x0); % Set 0 velocity first
params_new.psi = 1*(1-3*exp(-20.*params_new.sigma));
params_new.M = 5*10^-4/params_new.M0; 
params_new.N_terminus = rho_i*g*1/params_new.N0;
params_new.r = rho_i/rho_w;

QNS_new = [Q; N; S];
Qs_new = nan.*ones(params_new.Nt,params_new.Nx);
Ns_new = nan.*ones(params_new.Nt,params_new.Nx);
Ss_new = nan.*ones(params_new.Nt,params_new.Nx);

params_new.dt = endTime/(params_new.Nt*params_new.th0);
%% Solving
for t=2:params.Nt % same as paras_new
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

    eqns_new = @(QNS_new) hydro_eqns_new(QNS_new,params_new);
    [QNS_new,F,exitflag,output,JAC] = fsolve(eqns,QNS_new,options);
    Q_new = QNS_new(1:params.Nx);
    N_new = QNS_new(params.Nx+1:2*params.Nx);
    S_new = QNS_new(2*params.Nx+1:3*params.Nx);
    params_new.S_old_new = S_new;
    Qs_new(t,:)=Q_new;
    Ns_new(t,:)=N_new;
    Ss_new(t,:)=S_new;
end

%% plotting
Q_diff = abs(Qs_new-Qs);
N_diff = abs(Ns_new-Ns);
S_diff = abs(Ss_new-Ss);

ts = linspace(0,params.Nt*params.dt*params.th0./params.year,params.Nt);
subplot(3,1,1);surface(ts,params.sigma*params.x0,Q_diff',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Flow (m^3/s)');set(gca,'Ydir','Reverse')
subplot(3,1,2);surface(ts,params.sigma*params.x0,N_diff',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Effective Pressure (Pa)');set(gca,'Ydir','Reverse')
subplot(3,1,3);surface(ts,params.sigma*params.x0,S_diff',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Surface Area (m^2)');set(gca,'Ydir','Reverse')

figure('Name','Evolution of Drainage');
subplot(3,1,1);
plot(ts,N_diff(:,1),'DisplayName','Channel Entrance');
hold on;
plot(ts,N_diff(:,end),'--','DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)');
ylabel('Effective Pressure, \emph{N} (Pa)');
title('Effective Pressure over time');

subplot(3,1,2);
plot(ts,Q_diff(:,1),'DisplayName','Channel Entrance');
hold on;
plot(ts,Q_diff(:,end),'DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)');
ylabel('channel flow rate, \emph{Q} $(m^{3} s^{-1})$');
title('Flow over time');
legend('Location','northwest');

subplot(3,1,3);
plot(ts,S_diff(:,1),'-o','DisplayName','Channel Entrance');
hold on;
plot(ts,S_diff(:,end),'-o','DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)');
ylabel('channel cross-sectional area, \emph{S} $(m^2)$');
title('Channel area over time');
legend('Location','northwest');

%% Functions
function F = hydro_eqns(QNS,params)
    M = params.M;
    Q_in = 10/params.Q0;
    Nx = params.Nx;
    Q = QNS(1:Nx);
    N = QNS(Nx+1:2*Nx);
    S = QNS(2*Nx+1:3*Nx);
    % Q
    fq(1) = Q(1) - Q_in; % Boundary condition

    fq(2:Nx-1) = (params.eps_r*(params.r-1))*abs(Q(2:Nx-1)).^3./S(2:Nx-1).^(8/3) +...
                        params.eps_r.*S(2:Nx-1).*N(2:Nx-1).^3 + M + ...
                        params.u(2:Nx-1).*(S(3:Nx)-S(1:Nx-2))./(2*params.dsigma(2:Nx-1)) - ...
                        (Q(3:Nx)-Q(1:Nx-2))./(2*params.dsigma(2:Nx-1));
    fq(Nx) = (params.eps_r*(params.r-1))*abs(Q(Nx)).^3./S(Nx).^(8/3) +...
                        params.eps_r*S(Nx).*N(Nx).^3 + M + ...
                        params.u(Nx).*(S(Nx)-S(Nx-1))./(params.dsigma(Nx-1)) - ...
                        (Q(Nx)-Q(Nx-1))./(params.dsigma(Nx-1)); % one sided difference instead
    % N 
    fn(1) = Q(1).*abs(Q(1))./(S(1).^(8/3)) - params.psi(1) - ...
                        params.delta*(N(2)-N(1))./(params.dsigma(1)); % use 1 sided difference instead of symmetry argument
    fn(2:Nx-1) =  Q(2:Nx-1).*abs(Q(2:Nx-1))./(S(2:Nx-1).^(8/3)) - params.psi(2:Nx-1) - ...
                        params.delta*(N(3:Nx)-N(1:Nx-2))./(2*params.dsigma(2:Nx-1));
    fn(Nx) = N(Nx) - params.N_terminus; % Boundary condition

    % S
    fs(1) = abs(Q(1)).^3./(S(1).^(8/3)) - ... 
                        S(1).*N(1).^3 - ...
                        params.u(1).*(S(2)-S(1))./(params.dsigma(1)) - ...
                        (S(1)-params.S_old(1))./params.dt; % one sided differece
    fs(2:Nx-1)= abs(Q(2:Nx-1)).^3./(S(2:Nx-1).^(8/3)) - ... 
                        S(2:Nx-1).*N(2:Nx-1).^3 - ...
                        params.u(2:Nx-1).*(S(3:Nx)-S(1:Nx-2))./(2*params.dsigma(2:Nx-1)) - ...
                        (S(2:Nx-1)-params.S_old(2:Nx-1))./params.dt; 
    fs(Nx)= abs(Q(Nx)).^3./(S(Nx).^(8/3)) - ... 
                        S(Nx).*N(Nx).^3 - ...
                        params.u(Nx).*(S(Nx)-S(Nx-1))./(params.dsigma(Nx-1)) - ...
                        (S(Nx)-params.S_old(Nx))./params.dt; % one sided differece

    F = [fq;fn;fs];
end

function F = hydro_eqns_new(QNS_new,params)
    M = params.M;
    Q_in = 10/params.Q0;
    Nx = params.Nx;
    Q = QNS_new(1:Nx);
    N = QNS_new(Nx+1:2*Nx);
    S = QNS_new(2*Nx+1:3*Nx);
    % Q
    fq(1) = Q(1) - Q_in; % Boundary condition

    fq(2:Nx-1) = (params.eps_r*(params.r-1))*abs(Q(2:Nx-1)).^3./S(2:Nx-1).^(8/3) +...
                        params.eps_r.*S(2:Nx-1).*N(2:Nx-1).^3 + M + ...
                        params.beta.*params.u(2:Nx-1).*(S(3:Nx)-S(1:Nx-2))./(2*params.dsigma(2:Nx-1)) - ...
                        (Q(3:Nx)-Q(1:Nx-2))./(2*params.dsigma(2:Nx-1));
    fq(Nx) = (params.eps_r*(params.r-1))*abs(Q(Nx)).^3./S(Nx).^(8/3) +...
                        params.eps_r*S(Nx).*N(Nx).^3 + M + ...
                        params.beta.*params.u(Nx).*(S(Nx)-S(Nx-1))./(params.dsigma(Nx-1)) - ...
                        (Q(Nx)-Q(Nx-1))./(params.dsigma(Nx-1)); % one sided difference instead
    % N 
    fn(1) = Q(1).*abs(Q(1))./(S(1).^(8/3)) - params.psi(1) - ...
                        params.delta*(N(2)-N(1))./(params.dsigma(1)); % use 1 sided difference instead of symmetry argument
    fn(2:Nx-1) =  Q(2:Nx-1).*abs(Q(2:Nx-1))./(S(2:Nx-1).^(8/3)) - params.psi(2:Nx-1) - ...
                        params.delta*(N(3:Nx)-N(1:Nx-2))./(2*params.dsigma(2:Nx-1));
    fn(Nx) = N(Nx) - params.N_terminus; % Boundary condition

    % S
    fs(1) = abs(Q(1)).^3./(S(1).^(8/3)) - ... 
                        S(1).*N(1).^3 - ...
                        params.u(1).*(S(2)-S(1))./(params.dsigma(1)) - ...
                        (S(1)-params.S_old(1))./params.dt; % one sided differece
    fs(2:Nx-1)= abs(Q(2:Nx-1)).^3./(S(2:Nx-1).^(8/3)) - ... 
                        S(2:Nx-1).*N(2:Nx-1).^3 - ...
                        params.u(2:Nx-1).*(S(3:Nx)-S(1:Nx-2))./(2*params.dsigma(2:Nx-1)) - ...
                        (S(2:Nx-1)-params.S_old(2:Nx-1))./params.dt; 
    fs(Nx)= abs(Q(Nx)).^3./(S(Nx).^(8/3)) - ... 
                        S(Nx).*N(Nx).^3 - ...
                        params.u(Nx).*(S(Nx)-S(Nx-1))./(params.dsigma(Nx-1)) - ...
                        (S(Nx)-params.S_old(Nx))./params.dt; % one sided differece

    F = [fq;fn;fs];
end