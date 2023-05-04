clear all
close all
clc
%% Bed parameters
params.b0 = -100;           %bed topo at x=0
params.bx = -1e-3;          %linear bed slope

params.sill_min = 200e3;   %sill min x position
params.sill_max = 210e3;   %sill max x position
params.sill_slope = 1e-3;   %slope of sill

%% Physical parameters
params.year = 3600*24*365;  %number of seconds in a year
params.Aglen = 4.337e-25;   %ice softness parameter
params.nglen = 3;           %Glen's exponent
params.Bglen = params.Aglen^(-1/params.nglen);
params.m     = 1/params.nglen;  %sliding exponent (power law)
params.accum = 1/params.year; %SMB (constant here)
params.C     = 7e6;        %sliding coefficient (power law)
params.rhoi  = 917;         %ice density
params.rhow  = 1028;        %water density
params.g     = 9.81;        %gravity accel

%% Scaling params (coupled model equations solved in non-dim form)
params.hscale = 1000;               %thickness scaling
params.ascale = 0.1/params.year;    %SMB scaling 
params.uscale = (params.rhoi*params.g*params.hscale*params.ascale./params.C).^(1/(params.m+1)); %velocity scaling
params.xscale = params.uscale*params.hscale/params.ascale;  %horizontal distance scaling
params.tscale = params.xscale/params.uscale;                %time scaling
params.eps    = params.Bglen*((params.uscale/params.xscale)^(1/params.nglen))/(2*params.rhoi*params.g.*params.hscale);  %epsilon param (Schoof 2007)
params.lambda = 1 - (params.rhoi/params.rhow);  %density difference (lambda param Schoof 2007)
params.transient = 0;   %0 if solving for steady-state, 1 if solving for transient evolution

%% Grid parameters
params.tfinal = 10e3.*params.year;   %length of transient simulation
params.Nt = 1e2;                    %number of time steps
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
xg = 200e3/params.xscale;
hf = (-bed(xg.*params.xscale,params)./params.hscale)/(1-params.lambda);
h = 1 - (1-hf).*params.sigma;
u = 0.3*(params.sigma_elem.^(1/3)) + 1e-3;

params.h_old = h;
params.xg_old = xg;

sig_old = params.sigma;
sige_old = params.sigma_elem;
huxg0 = [h;u;xg];
options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3,'StepTolerance',1e-10,'FunctionTolerance',1e-10);
flf = @(huxg) flowline_eqns(huxg,params);

[huxg_init,F,exitflag,output,JAC] = fsolve(flf,huxg0,options);

h = huxg_init(1:params.Nx);
u = huxg_init(params.Nx+1:2*params.Nx);
xg = huxg_init(end);
hf = (-bed(xg.*params.xscale,params)/params.hscale)/(1-params.lambda);

u_bl = u_schoof(-bed(xg.*params.xscale,params)/(1-params.lambda),params);
u_num = u(end).*params.uscale;
num_err = abs((u_num-u_bl)/u_bl);
disp(['Error on initial S-S: ' num2str(100*num_err) '%']); %calculate departure of solution from Schoof 2007 BL approximation for GL flux)
save("SSA_simple_result.mat","h","u","xg","params");
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
figure(1);
ts = linspace(0,params.tfinal./params.year,params.Nt);
subplot(3,1,1);plot(ts,xgs.*params.xscale./1e3,'linewidth',3);xlabel('time (yr)');ylabel('x_g')
subplot(3,1,2);contourf(ts,params.sigma_elem,hs'.*params.hscale);colorbar;xlabel('time (yr)');xlabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse')
subplot(3,1,3);contourf(ts,params.sigma,us'.*params.uscale.*params.year);colorbar;xlabel('time (yr)');xlabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse')

%% Implicit system of equations function (using discretization scheme from Schoof 2007)
function F = flowline_eqns(huxg,params)
    
    %vars unpack
    h = huxg(1:params.Nx);
    u = huxg(params.Nx+1:2*params.Nx);
    xg = huxg(2*params.Nx+1);
    hf = (-bed(xg.*params.xscale,params)/params.hscale)/(1-params.lambda);

    %grid params unpack
    dt = params.dt/params.tscale;
    ds = params.dsigma;
    Nx = params.Nx;
    N1 = params.N1;
    sigma = params.sigma;
    sigma_elem = params.sigma_elem;
    b = -bed(xg.*sigma.*params.xscale,params)/params.hscale;
    Fh = zeros(Nx,1);
    Fu = zeros(Nx,1);
    
    %physical params unpack
    m     = params.m;
    nglen = params.nglen;
    lambda= params.lambda;
    accum = params.accum;
    a = accum/params.ascale;
    eps = params.eps;
    ss = params.transient;
       
    %previous time step unpack
    h_old = params.h_old;
    xg_old = params.xg_old;
    
    %thickness
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
    Fu(1)      = (4*eps).*(1./(xg.*ds(1)).^((1/nglen)+1)).*...
                 (h(2).*(u(2)-u(1)).*abs(u(2)-u(1)).^((1/nglen)-1) -...
                  h(1).*(2*u(1)).*abs(2*u(1)).^((1/nglen)-1)) -...
                  u(1).*abs(u(1)).^(m-1) -...
                  0.5.*(h(1)+h(2)).*(h(2)-b(2)-h(1)+b(1))./(xg.*ds(1));
    Fu(2:Nx-1) = (4*eps).*(1./(xg.*ds(2:Nx-1)).^((1/nglen)+1)).*...
                 (h(3:Nx).*(u(3:Nx)-u(2:Nx-1)).*abs(u(3:Nx)-u(2:Nx-1)).^((1/nglen)-1) -...
                  h(2:Nx-1).*(u(2:Nx-1)-u(1:Nx-2)).*abs(u(2:Nx-1)-u(1:Nx-2)).^((1/nglen)-1)) -...
                  u(2:Nx-1).*abs(u(2:Nx-1)).^(m-1) -...
                  0.5.*(h(2:Nx-1)+h(3:Nx)).*(h(3:Nx)-b(3:Nx)-h(2:Nx-1)+b(2:Nx-1))./(xg.*ds(2:Nx-1));
    Fu(N1)     = (u(N1+1)-u(N1))/ds(N1) - (u(N1)-u(N1-1))/ds(N1-1);
    Fu(Nx)     = (1./(xg.*ds(Nx-1)).^(1/nglen)).*...
                 (abs(u(Nx)-u(Nx-1)).^((1/nglen)-1)).*(u(Nx)-u(Nx-1)) - lambda*hf/(8*eps);
             
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