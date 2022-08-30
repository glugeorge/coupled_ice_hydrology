clear all
close all
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

params.psi0 = params.rho_i*params.g*params.h0/params.x0;
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
params.tfinal = 10000.*params.year;   %length of transient simulation
params.Nt = 100;                    %number of time steps
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

%% Solve for steady-state initial conditions
params.accum = 1/params.year;
xg = 10e3/params.x0;
hf = (-bed(xg.*params.x0,params)/params.h0)/params.r;
h = 1 - (1-hf).*params.sigma;
u = 0.3*(params.sigma_elem.^(1/3)) + 1e-3;

params.h_old = h;
params.xg_old = xg;

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


%% Calculate transient GL evolution over bedrock peak
xgs = nan.*ones(1,params.Nt);
hs = nan.*ones(params.Nt,params.Nx);
us = nan.*ones(params.Nt,params.Nx);
huxg_t = huxg_init;
params.h_old = huxg_t(1:params.Nx);
params.xg_old = huxg_t(end);

params.transient = 1;
params.accum = 0.8/params.year;
for t=1:params.Nt
    flf = @(huxg) flowline_eqns(huxg,params);
    [huxg_t,F,exitflag,output,JAC] = fsolve(flf,huxg_t,options);
    
    params.h_old = huxg_t(1:params.Nx);
    params.xg_old = huxg_t(end);
    t
    xgs(t) = huxg_t(end);
    hs(t,:) = huxg_t(1:params.Nx)';
    us(t,:) = huxg_t(params.Nx+1:end-1)';
end


%% Plot transient solution
figure(2);
ts = linspace(0,params.tfinal./params.year,params.Nt);
subplot(3,1,1);plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time (yr)');ylabel('x_g')
%subplot(3,1,2);contourf(ts,params.sigma_elem,hs'.*params.hscale);colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse')
subplot(3,1,3);contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse')

subplot(3,1,2);surface(ts,params.sigma_elem,hs'.*params.h0,EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse')
%subplot(3,1,3);surface(ts,params.sigma,us'.*params.uscale.*params.year,EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse')
%%
N = 100000;
u_new = us'.*params.u0;
tau = params.C.*N.*(u_new./(u_new+params.As*(params.C*N)^params.n)); 
figure(3)
contourf(ts,params.sigma,tau);colorbar;xlabel('time (yr)');ylabel('sigma');title('shear stress (Pa)');set(gca,'Ydir','Reverse');
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
    N = 100000*ones(size(h_old))./params.N0;
    
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

%% Schoof GL function

function us = u_schoof(hg,params)

    us = (((params.Aglen*(params.rhoi*params.g)^(params.nglen+1) * params.lambda^params.nglen)/(4^params.nglen * params.C))^(1/(params.m+1)))*(hg)^(((params.m+params.nglen+3)/(params.m+1))-1);
    
end
