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
params.A = 4.227e-25; % From Alex Robel's code
params.n = 3;
params.rho_i = 917;
params.rho_w = 1028;
params.g = 9.81;
params.C = 0.84*0.5; % From Hewitt and Fowler, 2008
params.As = 2*2.4e-24/(0.5*params.C^params.n); % Calculated 
params.p = 1/3; % From Hewitt's Karthaus slides 
params.q = 1/3; % From Hewitt's Karthaus slides
params.f = 0.07; % From Kingslake thesis
params.K0 = 10^-24; % From Kingslake thesis  
params.L = 3.3e5; % Kingslake thesis
params.year = 3600*24*365;
%% Scaling params (coupled model equations solved in non-dim form)
params.x0 = 10*10^3;
params.h0 = 100;
params.Q0 = 1500;

params.psi0 = params.rho_w*params.g*abs(params.bx);
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
params.Nx = 200;                    %number of grid points
params.N1 = 100;                    %number of grid points in coarse domain
params.sigGZ = 0.97;                %extent of coarse grid (where GL is at sigma=1)
sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
params.sigma = [sigma1, sigma2(2:end)]';    %grid points on velocity (includes GL, not ice divide)
params.sigma_elem = [0;(params.sigma(1:params.Nx-1) + params.sigma(2:params.Nx))./2]; %grid points on thickness (includes divide, not GL)
params.dsigma = diff(params.sigma); %grid spacing

%% Grid parameters - hydro
params.sigma_h = linspace(0,1,params.Nx)';
params.dsigma_h = diff(params.sigma_h); %grid spacing

%% Establish timings
params.year = 3600*24*365;  %number of seconds in a year
params.Nt = 100;                    %number of time steps
params.end_year = 100;

params.dt = params.end_year*params.year/params.Nt;

%% Initial "steady state" ice conditions
params.accum = 1/params.year;
xg = 10e3/params.x0;
hf = (-bed(xg.*params.x0,params)/params.h0)/params.r;
h = 1 - (1-hf).*params.sigma;
u = 0.3*(params.sigma_elem.^(1/3)) + 1e-3;

params.h_old = h;
params.xg_old = xg;
params.N_ice = 100000*ones(size(params.h_old))./params.N0;
sig_old = params.sigma;
sige_old = params.sigma_elem;
huxg0 = [h;u;xg];
options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
flf = @(huxg) flowline_eqns(huxg,params);

[huxg_init,F,exitflag,output,JAC] = fsolve(flf,huxg0,options);

h = huxg_init(1:params.Nx);
u = huxg_init(params.Nx+1:2*params.Nx);
xg = huxg_init(end);
hf = (-bed(xg.*params.x0,params)/params.h0)/(params.r);

%% Now for evolution 
Q = 0.001*ones(params.Nx,1);
N = ones(params.Nx,1);
S = 5/params.S0*ones(params.Nx,1); 
params.S_old = S;
params.M = 5*10^-4/params.M0; 
params.N_terminus = 0;%rho_i*g*h_interp(end)*params.h0/params.N0;

QNS = [Q; N; S];

Qs = nan.*ones(params.Nt,params.Nx);
Ns = nan.*ones(params.Nt,params.Nx);
Ss = nan.*ones(params.Nt,params.Nx);
xgs = nan.*ones(1,params.Nt);
hs = nan.*ones(params.Nt,params.Nx);
us = nan.*ones(params.Nt,params.Nx);
huxg_t = huxg_init;
params.h_old = huxg_t(1:params.Nx);
params.xg_old = huxg_t(end);
xgs(1) = params.xg_old;
params.transient = 1;
params.accum = 0.8/params.year;

params.h_old = h;
params.xg_old = xg;
params.S_old = S;
params.transient = 1;
params.accum = 0.8/params.year;
params.transient_h = 0; % 0 means use steady state

for t=2:params.Nt
    t
    % Update ice
    params.xg_old = xgs(t-1);

    flf = @(huxg) flowline_eqns(huxg,params);
    [huxg_t,F,exitflag,output,JAC] = fsolve(flf,huxg_t,options);
    params.h_old = huxg_t(1:params.Nx);
    xgs(t) = huxg_t(end);
    params.xg = xgs(t);
    hs(t,:) = huxg_t(1:params.Nx)';
    us(t,:) = huxg_t(params.Nx+1:end-1)';

    % Update params from ice that are used in hydro
    h_grid = params.sigma_elem*xgs(t); % nondimensionalized xgrid for h
    u_grid = params.sigma*xgs(t); % nondimensionalized xgrid for u
    h_interp = interp1(h_grid,hs(t,:),params.sigma_h*params.xg,'linear','extrap');
    u_interp = interp1(u_grid,us(t,:),params.sigma_h*params.xg,'linear','extrap');
    params.u = u_interp;
    % define psi with realistic h
    phi_b = 0.001; % slope
    p = params.rho_i*params.g*h_interp*params.h0;
    params.psi = (params.rho_w*params.g*sin(phi_b)-gradient(p)./gradient(params.sigma_h.*params.x0.*params.xg))./params.psi0;
    
    % Solve for hydro using previous step's ice
    
    eqns = @(QNS) hydro_eqns(QNS,params);
    [QNS,F,exitflag,output,JAC] = fsolve(eqns,QNS,options);
    Q = QNS(1:params.Nx);
    N = QNS(params.Nx+1:2*params.Nx);
    S = QNS(2*params.Nx+1:3*params.Nx);
    params.S_old = S;
    Qs(t,:)=Q;
    Ns(t,:)=N;
    Ss(t,:)=S;
    params.N_ice = interp1(params.sigma_h,N,params.sigma_elem,'linear','extrap');
end

%% Plotting
ts = linspace(0,params.end_year,params.Nt);
figure();
subplot(3,1,1);plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time (yr)');ylabel('x_g');
subplot(3,1,2);surface(ts,params.sigma_elem,hs'.*params.h0,EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse');
subplot(3,1,3);contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse');

figure();
Q_dim = Qs.*params.Q0;
N_dim = Ns.*params.N0;
S_dim = Ss.*params.S0;
subplot(3,1,1);surface(ts,params.sigma_h,Q_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Flow (m^3/s)');set(gca,'Ydir','Reverse')
subplot(3,1,2);surface(ts,params.sigma_h,N_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Effective Pressure (Pa)');set(gca,'Ydir','Reverse')
subplot(3,1,3);surface(ts,params.sigma_h,S_dim',EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('distance');title('Surface Area (m^2)');set(gca,'Ydir','Reverse')

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
%% Functions
function F = flowline_eqns(huxg,params)
    
    %vars unpack
    h = huxg(1:params.Nx);
    u = huxg(params.Nx+1:2*params.Nx);
    xg = huxg(2*params.Nx+1);
    hf = (-bed(xg.*params.x0,params)/params.h0)/(params.r);

    %grid params unpack
    dt = params.dt/params.t0;
    ds = params.dsigma;
    Nx = params.Nx;
    N1 = params.N1;
    sigma = params.sigma;
    sigma_elem = params.sigma_elem;
    b = -bed(xg.*sigma.*params.x0,params)/params.h0;
    Fh = zeros(Nx,1);
    Fu = zeros(Nx,1);
    
    %physical params unpack
    m     = 1/params.n;
    nglen = params.n;
    accum = params.accum;
    a = accum/params.a0;
    gamma = params.gamma;
    alpha = params.alpha;
    ss = params.transient;
       
    %previous time step unpack
    h_old = params.h_old;
    xg_old = params.xg_old;
    N = params.N_ice;
    
    %thickness - stays same
    Fh(1)      = ss.*(h(1)-h_old(1))./dt + (2.*h(1).*u(1))./(ds(1).*xg) - a;
    Fh(2)      = ss.*(h(2)-h_old(2))./dt -...
                    ss.*sigma_elem(2).*(xg-xg_old).*(h(3)-h(1))./(2*dt.*ds(2).*xg) +...
                        (h(2).*(u(2)+u(1)))./(2*xg.*ds(2)) -...
                            a;
    Fh(3:Nx-1) = ss.*(h(3:Nx-1)-h_old(3:Nx-1))./dt -...
                    ss.*sigma_elem(3:Nx-1).*(xg-xg_old).*(h(4:Nx)-h(2:Nx-2))./(2*dt.*ds(3:Nx-1).*xg) +...
                        (h(3:Nx-1).*(u(3:Nx-1)+u(2:Nx-2)) - h(2:Nx-2).*(u(2:Nx-2)+u(1:Nx-3)))./(2*xg.*ds(3:Nx-1)) -...
                            a;
    
    Fh(N1) = (1+0.5*(1+(ds(N1)/ds(N1-1))))*h(N1) - 0.5*(1+(ds(N1)/ds(N1-1)))*h(N1-1) - h(N1+1);
                        
    Fh(Nx)     = ss.*(h(Nx)-h_old(Nx))./dt -...
                    ss.*sigma_elem(Nx).*(xg-xg_old).*(h(Nx)-h(Nx-1))./(dt.*ds(Nx-1).*xg) +...
                        (h(Nx).*(u(Nx)+u(Nx-1)) - h(Nx-1).*(u(Nx-1)+u(Nx-2)))./(2*xg.*ds(Nx-1)) -...
                            a;

	%velocity
    Fu(1)      = (alpha).*(1./(xg.*ds(1)).^((1/nglen)+1)).*...
                 (h(2).*(u(2)-u(1)).*abs(u(2)-u(1)).^((1/nglen)-1) -...
                  h(1).*(2*u(1)).*abs(2*u(1)).^((1/nglen)-1)) -...
                  gamma.*N(1).* u(1)./(u(1) + N(1).^nglen) -...
                  0.5.*(h(1)+h(2)).*(h(2)-b(2)-h(1)+b(1))./(xg.*ds(1));
    Fu(2:Nx-1) = (alpha).*(1./(xg.*ds(2:Nx-1)).^((1/nglen)+1)).*...
                 (h(3:Nx).*(u(3:Nx)-u(2:Nx-1)).*abs(u(3:Nx)-u(2:Nx-1)).^((1/nglen)-1) -...
                  h(2:Nx-1).*(u(2:Nx-1)-u(1:Nx-2)).*abs(u(2:Nx-1)-u(1:Nx-2)).^((1/nglen)-1)) -...
                  gamma.*N(2:Nx-1).*u(2:Nx-1)./(u(2:Nx-1) + N(2:Nx-1).^nglen) -...
                  0.5.*(h(2:Nx-1)+h(3:Nx)).*(h(3:Nx)-b(3:Nx)-h(2:Nx-1)+b(2:Nx-1))./(xg.*ds(2:Nx-1));
    Fu(N1)     = (u(N1+1)-u(N1))/ds(N1) - (u(N1)-u(N1-1))/ds(N1-1);
    Fu(Nx)     = alpha.*(1./(xg.*ds(Nx-1)).^(1/nglen)).*...
                 (abs(u(Nx)-u(Nx-1)).^((1/nglen)-1)).*(u(Nx)-u(Nx-1)) - 0.5*(1-params.r)*hf;
             
    Fxg        = 3*h(Nx) - h(Nx-1) - 2*hf;         
    
    F = [Fh;Fu;Fxg];
end

function F = hydro_eqns(QNS,params)
    M = params.M;
    Q_in = 10/params.Q0;
    Nx = params.Nx;
    Q = QNS(1:Nx);
    N = QNS(Nx+1:2*Nx);
    S = QNS(2*Nx+1:3*Nx);
    xg = params.xg;
    xg_old = xg;
    sigma = params.sigma_h;
    dt = params.dt/params.th0;
    ssh = params.transient_h;

    % Q
    fq(1) = Q(1) - Q_in; % Boundary condition

    fq(2:Nx-1) = (params.eps_r*(params.r-1))*abs(Q(2:Nx-1)).^3./S(2:Nx-1).^(8/3) +...
                        params.eps_r.*S(2:Nx-1).*N(2:Nx-1).^3 + M + ...
                        params.eps_r.*params.beta.*params.u(2:Nx-1).*(S(3:Nx)-S(1:Nx-2))./(2*xg*params.dsigma_h(2:Nx-1)) - ...
                        (Q(3:Nx)-Q(1:Nx-2))./(2*xg*params.dsigma_h(2:Nx-1));
    fq(Nx) = (params.eps_r*(params.r-1))*abs(Q(Nx)).^3./S(Nx).^(8/3) +...
                        params.eps_r*S(Nx).*N(Nx).^3 + M + ...
                        params.eps_r.*params.beta.*params.u(Nx).*(S(Nx)-S(Nx-1))./(xg*params.dsigma_h(Nx-1)) - ...
                        (Q(Nx)-Q(Nx-1))./(xg*params.dsigma_h(Nx-1)); % one sided difference instead
    % N 
    fn(1) = Q(1).*abs(Q(1))./(S(1).^(8/3)) - params.psi(1) - ...
                        params.delta*(N(2)-N(1))./(xg*params.dsigma_h(1)); % use 1 sided difference instead of symmetry argument
    fn(2:Nx-1) =  Q(2:Nx-1).*abs(Q(2:Nx-1))./(S(2:Nx-1).^(8/3)) - params.psi(2:Nx-1) - ...
                        params.delta*(N(3:Nx)-N(1:Nx-2))./(2*xg*params.dsigma_h(2:Nx-1));
    fn(Nx) = N(Nx) - params.N_terminus; % Boundary condition

    % S
    fs(1) = abs(Q(1)).^3./(S(1).^(8/3)) - ... 
                        S(1).*N(1).^3 + ...
                        (ssh*sigma(1)*(xg-xg_old)/dt - params.beta.*params.u(1)).*(S(2)-S(1))./(xg*params.dsigma_h(1)) - ...
                        ssh*(S(1)-params.S_old(1))./dt; % one sided differece
    fs(2:Nx-1)= abs(Q(2:Nx-1)).^3./(S(2:Nx-1).^(8/3)) - ... 
                        S(2:Nx-1).*N(2:Nx-1).^3 + ...
                        (ssh*sigma(2:Nx-1).*(xg-xg_old)./dt - params.beta.*params.u(2:Nx-1)).*(S(3:Nx)-S(1:Nx-2))./(2*xg*params.dsigma_h(2:Nx-1)) - ...
                        ssh*(S(2:Nx-1)-params.S_old(2:Nx-1))./dt; 
    fs(Nx)= abs(Q(Nx)).^3./(S(Nx).^(8/3)) - ... 
                        S(Nx).*N(Nx).^3 + ...
                        (ssh*sigma(Nx).*(xg-xg_old)./dt - params.beta.*params.u(Nx)).*(S(Nx)-S(Nx-1))./(xg*params.dsigma_h(Nx-1)) - ...
                        ssh*(S(Nx)-params.S_old(Nx))./dt; % one sided differece

    F = [fq;fn;fs];
end

%% Bed topography function
function b = bed(x,params)
    xsill = x>params.sill_min & x<params.sill_max;
    xdsill = x>=params.sill_max;
    sill_length = params.sill_max-params.sill_min;
    
    b = params.b0 + params.bx.*x;
    
    b(xsill) = params.b0 + (params.bx*params.sill_min) + params.sill_slope.*(x(xsill)-params.sill_min);

    b(xdsill) = params.b0 + (params.bx*params.sill_min) + params.sill_slope.*sill_length + ...
            params.bx*(x(xdsill)-params.sill_max);

end

