clear all
close all
clc
%% For using same initial condition as coupled 

filename = 'C_0.1_A_4.9_c.mat';
load(filename);
% case where N is the same as coupled, kept constant
N = results.init_cond(results.params.Nx+1:2*results.params.Nx);
h = results.init_cond(3*results.params.Nx+1:4*results.params.Nx);
u = results.init_cond(4*results.params.Nx+1:5*results.params.Nx);
xg = results.init_cond(5*results.params.Nx+1);
params.sigma_h = results.params.sigma_h;
params.fixed_N_grid = params.sigma_h*xg; 
params.N_scaled = N*results.params.N0;
h_scaled = h*results.params.h0;
u_scaled = u*results.params.u0;

% case where N is unifrom and constant and C varies
% params.N_scaled = ones(size(N)).*100000;
% load C_new.mat
% params.C_new = C_new;
%% Bed parameters
params.b0 = -100;           %bed topo at x=0
params.bx = -1e-3;          %linear bed slope

params.sill_min = 2000e3;   %sill min x position
params.sill_max = 2100e3;   %sill max x position
params.sill_slope = 1e-3;   %slope of sill

%% Physical parameters
params.A = 2.9e-25; % From Alex Robel's code
params.n = 3;
params.rho_i = 917;
params.rho_w = 1028;
params.g = 9.81;
params.C = results.params.C; % 
params.As = 2.26e-21; % Calculated 
params.f = 0.07; % From Kingslake thesis
params.K0 = 10^-24; % From Kingslake thesis  
params.L = 3.3e5; % Kingslake thesis
params.year = 3600*24*365;
%% Scaling params (coupled model equations solved in non-dim form)
params.x0 = 10*10^3;
params.h0 = 100;
params.Q0 = 1500;

params.psi0 = params.rho_i*params.g*params.h0/params.x0/10;
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

params.transient = 0;   %0 if solving for steady-state, 1 if solving for transient evolution

%% Grid parameters
params.tfinal = 7.5e2.*params.year;   %length of transient simulation
params.Nt = 150;                    %number of time steps
params.dt = params.tfinal/params.Nt;%time step length
params.Nx = 200;                    %number of grid points
params.N1 = 100;                    %number of grid points in coarse domain
params.sigGZ = 0.97;                %extent of coarse grid (where GL is at sigma=1)
% params.sigma = linspace(0,1-(1/(2*params.Nx)),params.Nx)';
% params.sigma = [linspace(0,0.97,params.Nx/2)';linspace(0.97+(0.03/(params.Nx/2)),1-(0.03/(2*params.Nx)),params.Nx/2)']; %piecewise refined grid, with 30x finer resolution near GL
sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
params.sigma = [sigma1, sigma2(2:end)]';    %grid points on velocity (includes GL, not ice divide)
params.sigma_elem = [0;(params.sigma(1:params.Nx-1) + params.sigma(2:params.Nx))./2]; %grid points on thickness (includes divide, not GL)
params.dsigma = diff(params.sigma); %grid spacing


%% Test to make sure flowline equations make sense
params.accum = results.params.accum;%1/params.year;
%hf = (-bed(xg.*params.x0,params)/params.h0)/params.r;
h = h_scaled/params.h0;%1 - (1-hf).*params.sigma;
u = u_scaled/params.u0;%0.3*(params.sigma_elem.^(1/3)) + 1e-3;
params.h_old = h;
params.xg_old = xg;

sig_old = params.sigma;
sige_old = params.sigma_elem;
huxg0 = [h;u;xg];

res = flowline_eqns(huxg0,params);
h_diff = res(1:params.Nx);
u_diff = res(params.Nx+1:2*params.Nx);
xg_diff = res(end);

%% Solve for steady-state initial conditions

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
flf = @(huxg) flowline_eqns(huxg,params);

[huxg_init,F,exitflag,output,JAC] = fsolve(flf,huxg0,options);

h = huxg_init(1:params.Nx);
u = huxg_init(params.Nx+1:2*params.Nx);
xg = huxg_init(end);
hf = (-bed(xg.*params.x0,params)/params.h0)/(params.r);

%% Calculate steady state solution
%params.accum = 1/params.year;
params.A = results.params.A; 
params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));
flf = @(huxg) flowline_eqns(huxg,params);
[huxg_final,F,exitflag,output,JAC] = fsolve(flf,huxg0,options);
xg_f = huxg_final(end);

%% Calculate transient GL evolution over bedrock peak
xgs = nan.*ones(1,params.Nt);
hs = nan.*ones(params.Nt,params.Nx);
us = nan.*ones(params.Nt,params.Nx);
huxg_t = huxg_init;
params.h_old = huxg_t(1:params.Nx);
params.xg_old = huxg_t(end);
hs(1,:) = huxg_t(1:params.Nx);
us(1,:) = huxg_t(params.Nx+1:2*params.Nx);
xgs(1) = huxg_t(end);
params.transient = 1;
time_to_ss = 0;
for t=2:params.Nt
    flf = @(huxg) flowline_eqns(huxg,params);
    [huxg_t,F,exitflag,output,JAC] = fsolve(flf,huxg_t,options);
    
    params.h_old = huxg_t(1:params.Nx);
    params.xg_old = huxg_t(end);
    t
    xgs(t) = huxg_t(end);
    hs(t,:) = huxg_t(1:params.Nx)';
    us(t,:) = huxg_t(params.Nx+1:end-1)';
    if abs(xg_f - xgs(t)) < 0.001*xg_f && time_to_ss == 0
        time_to_ss = (t-1)*params.dt/params.year;
    end
end


%% Plot transient solution
figure();
ts = linspace(0,params.tfinal./params.year,params.Nt);
subplot(3,1,1);plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time (yr)');ylabel('x_g');hold on;plot(results.ts,results.xgs.*params.x0./1e3,'linewidth',3)
%subplot(3,1,2);contourf(ts,params.sigma_elem,hs'.*params.hscale);colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse')
subplot(3,1,3);contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse')

subplot(3,1,2);contourf(ts,params.sigma_elem,hs'.*params.h0);colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse')
%subplot(3,1,3);surface(ts,params.sigma,us'.*params.uscale.*params.year,EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse')
%% Plot evolution of N
figure();
hold on;
colors = ['r','y','g','b','k'];
for i=1:5
    plot(params.sigma*xgs(i)*params.x0/1000,interp1(params.fixed_N_grid,params.N_scaled,params.sigma*xgs(i)),"Color",colors(i))
    plot(params.sigma(end)*xgs(i)*params.x0/1000,interp1(params.fixed_N_grid,params.N_scaled,params.sigma(end)*xgs(i)),'o',"Color",colors(i))

end

%% Save results
results_constN.params = params;
results_constN.init_cond = huxg_init;
results_constN.steady_state = huxg_final;
results_constN.xgs = xgs;
results_constN.ts = ts;
results_constN.hs = hs';
results_constN.us = us';
results_constN.time_to_ss = time_to_ss; 
fname = strrep(filename,'c','uc');
%save(fname,'results_constN');

%% Implicit system of equations function (using discretization scheme from Schoof 2007)
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
    N = interp1(params.fixed_N_grid,params.N_scaled./params.N0,sigma*xg);
    
    % If using variable C
    %C_adjust = params.C_new./params.C;
    % else
    C_adjust = ones(size(N));

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
                  gamma.*C_adjust(1).*N(1).* (u(1)./(u(1) + (C_adjust(1).*N(1)).^nglen)).^m -...
                  0.5.*(h(1)+h(2)).*(h(2)-b(2)-h(1)+b(1))./(xg.*ds(1));
    Fu(2:Nx-1) = (alpha).*(1./(xg.*ds(2:Nx-1)).^((1/nglen)+1)).*...
                 (h(3:Nx).*(u(3:Nx)-u(2:Nx-1)).*abs(u(3:Nx)-u(2:Nx-1)).^((1/nglen)-1) -...
                  h(2:Nx-1).*(u(2:Nx-1)-u(1:Nx-2)).*abs(u(2:Nx-1)-u(1:Nx-2)).^((1/nglen)-1)) -...
                  gamma.*C_adjust(2:Nx-1).*N(2:Nx-1).*(u(2:Nx-1)./(u(2:Nx-1) + (C_adjust(2:Nx-1).*N(2:Nx-1)).^nglen)).^m -...
                  0.5.*(h(2:Nx-1)+h(3:Nx)).*(h(3:Nx)-b(3:Nx)-h(2:Nx-1)+b(2:Nx-1))./(xg.*ds(2:Nx-1));
    Fu(N1)     = (u(N1+1)-u(N1))/ds(N1) - (u(N1)-u(N1-1))/ds(N1-1);
    Fu(Nx)     = alpha.*(1./(xg.*ds(Nx-1)).^(1/nglen)).*...
                 (abs(u(Nx)-u(Nx-1)).^((1/nglen)-1)).*(u(Nx)-u(Nx-1)) - 0.5*(1-params.r)*hf;
             
    Fxg        = 3*h(Nx) - h(Nx-1) - 2*hf;         
    
    F = [Fh;Fu;Fxg];
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
