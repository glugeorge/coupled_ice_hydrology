S_initial = 10; % initial channel size [m^2]
Qin = 1;        % discharge into the top end of the channel  [m^3/s]  
dt = 0.01;      % time step [non-dim time]

Sdiff = differenceAfterOneTimeStep(0.01,1,1);

plot(Sdiff)
xlabel 'grid point number'
ylabel 'fractional difference between S solutions'


function[Sdiff] = differenceAfterOneTimeStep(dt,Qin,S_initial)

%% 1. Setup
params.Nt = 3;                    %number of time steps
params.year = 3600*24*365;        %number of seconds in a year
params.dt = dt;%;0.001;
params.Nx = 200;                  % number of grid points

T = params.Nt*params.dt;          % total time

t = linspace(0,T,params.Nt);      % time vector


params.sigma = linspace(0,1,params.Nx)';%[linspace(params.sigGZ/(params.Nx+0.5),params.sigGZ,params.Nx)]';
params.dsigma = diff(params.sigma); %grid spacing

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
% provide a initial guess of  N and Q
Q = 0.001*ones(params.Nx,1);
N = ones(params.Nx,1);
S = S_initial/params.S0*ones(params.Nx,1); 
params.u = zeros(params.Nx,1);
params.psi = params.psi*(1-3*exp(-20.*params.sigma));

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);

params.Qin =Qin/params.Q0;


%% 2. compute the new N, S and Q from the initial S using the implicit method
QNS = [Q; N; S];
params.S_old = S;
eqns = @(QNS) hydro_eqns(QNS,params);
[QNS,F,exitflag,output,JAC] = fsolve(eqns,QNS,options);
Q_imp = QNS(1:params.Nx);
N_imp = QNS(params.Nx+1:2*params.Nx);
S_imp = QNS(2*params.Nx+1:3*params.Nx);



%% 3. compute the new N and Q from the BVP method and use it to update S

% compute initial N and Q
opts = bvpset('RelTol',0.00000001,'AbsTol',0.0000001,'Stats','on');
solinit= bvpinit(params.sigma, @(x) Nyeguess(x,params));
dSdx = gradient(S(:,1))./gradient(params.sigma)'; % using previous time step's dSdx in coupled solver
sol = bvp5c(@(x,y) Nye_NQ_simple(x,y, params.sigma, S, params),...
    @(ya,yb) bc_N(ya,yb,params.Qin,params.N_terminus),...
    solinit, opts);
N = interp1(sol.x,sol.y(1,:),params.sigma);
Q = interp1(sol.x,sol.y(2,:),params.sigma);

% use N and Q to step forward S
S_bvp = S + params.dt.*(abs(Q.^3)./S.^(8/3) - S.*N.^3 - params.u(:,1).*dSdx);

% % find new N and Q
% dSdx = gradient(S(:,1))./gradient(params.sigma); % using previous time step's dSdx in coupled solver
% sol = bvp5c(@(x,y) Nye_NQ(x,y, params.sigma, S_bvp, params.u(:,1),dSdx, P),...
%     @(ya,yb) bc_N(ya,yb,params.Qin,params.N_terminus),...
%     solinit, opts);
% N_bvp = interp1(sol.x,sol.y(1,:),params.sigma)';
% Q_bvp = interp1(sol.x,sol.y(2,:),params.sigma)';


%% 4. compare
Sdiff = S_imp-S_bvp;

end


function F = hydro_eqns(QNS,params)
    M = params.M;
    Q_in = params.Qin; 
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

function dydx = Nye_NQ_simple(x,y,xmesh,S_i,params) % equation to solve


S_x = interp1(xmesh,S_i,x);
psi_x = interp1(xmesh,params.psi,x);
% u_x = interp1(xmesh,u,x);
% dSdx = interp1(xmesh,dSdx,x);
S83 = S_x.^(8/3);
N = y(1,:);
Q = y(2,:);
dydx = [(Q.*abs(Q)./S83 - psi_x )/params.delta % change between constant/variable psi
    params.eps_r*(params.r-1)*abs(Q).^3./S83 + params.eps_r*S_x.*N.^3 + params.M ]; % Use just P.M for efficiency
end

% Function for boundary conditions
function res = bc_N(ya,yb,Q_in,Nt) 
res = [ya(2)-Q_in
    yb(1)-Nt];
end
