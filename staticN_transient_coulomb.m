clear all
close all
clc
%% For using same initial condition as coupled 

load coulomb_retreat_5000yr.mat
params = results.params;
% case where N is the same as coupled, kept constant
N = results.init_cond(params.Nh+1:2*params.Nh);
h = results.init_cond(params.ice_start+1:params.ice_start+ params.Nx);
u = results.init_cond(params.ice_start + params.Nx+1:params.ice_start+2*params.Nx);
xg = results.init_cond(params.ice_start+2*params.Nx+1);

params.A = 0.9e-25;
params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));
params.fixed_N_grid = params.sigma_h*xg; 
params.N_scaled = N*params.N0;
h_scaled = h*results.params.h0;
u_scaled = u*results.params.u0;
params.transient = 0;

%% Test to make sure flowline equations make sense
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
hf = (-bed_schoof(xg.*params.x0,params)/params.h0)/(params.r);

%% Calculate steady state solution
%params.accum = 1/params.year;
params.A = 4.9e-25; 
params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));
flf = @(huxg) flowline_eqns(huxg,params);
[huxg_final,F,exitflag,output,JAC] = fsolve(flf,huxg_init,options);
xg_f = huxg_final(end);


%% Establish timings
params.year = 3600*24*365;  %number of seconds in a year
params.Nt =100;                    %number of time steps - normally 150
params.end_year = 5000; %normally 7500

params.dt = params.end_year*params.year/params.Nt;

%% Calculate transient GL evolution over bedrock peak
huxg_t = huxg0;
xgs = nan.*ones(1,params.Nt);
hs = nan.*ones(params.Nt,params.Nx);
us = nan.*ones(params.Nt,params.Nx);
params.h_old = huxg_t(1:params.Nx);
params.xg_old = xg;
hs(1,:) = huxg_t(1:params.Nx);
us(1,:) = huxg_t(params.Nx+1:2*params.Nx);
xgs(1) = huxg_t(end);
params.transient = 1;
time_to_ss = 0;
for t=1:params.Nt
    flf = @(huxg) flowline_eqns(huxg,params);
    [huxg_t,F,exitflag,output,JAC] = fsolve(flf,huxg_t,options);
    
    params.h_old = huxg_t(1:params.Nx);
    params.xg_old = huxg_t(end);
    t
    xgs(t) = huxg_t(end);
    hs(t,:) = huxg_t(1:params.Nx)';
    us(t,:) = huxg_t(params.Nx+1:end-1)';
    %if abs(xg_f - xgs(t)) < 0.001*xg_f && time_to_ss == 0
    %    time_to_ss = (t-1)*params.dt/params.year;
    %end
end


%% Plot transient solution
figure();
ts = linspace(0,params.end_year,params.Nt);
subplot(3,1,1);plot(ts,xgs.*params.x0./1e3,'linewidth',3);xlabel('time (yr)');ylabel('x_g');hold on;plot(results.ts,results.xgs.*params.x0./1e3,'linewidth',3)
%subplot(3,1,2);contourf(ts,params.sigma_elem,hs'.*params.hscale);colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse')
subplot(3,1,3);contourf(ts,params.sigma,us'.*params.u0.*params.year);colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse')

subplot(3,1,2);contourf(ts,params.sigma_elem,hs'.*params.h0);colorbar;xlabel('time (yr)');ylabel('sigma');title('thickness (m)');set(gca,'Ydir','Reverse')
%subplot(3,1,3);surface(ts,params.sigma,us'.*params.uscale.*params.year,EdgeColor='None');colorbar;xlabel('time (yr)');ylabel('sigma');title('velocity (m/yr)');set(gca,'Ydir','Reverse')
%% Plot evolution of N
figure();
pointsize = 10;

scatter(params.sigma(end).*xgs.*params.x0/1000,interp1(params.fixed_N_grid,params.N_scaled,params.sigma(end)*xgs),pointsize,ts); a = colorbar();colormap jet; ylabel(a,'time (yr)')
hold on;
plot(params.fixed_N_grid*params.x0/1000,params.N_scaled);ylabel('N (Pa)');xlabel('distance from divide (km)');

%% Save results
results_constN.params = params;
results_constN.init_cond = huxg_init;
%results_constN.steady_state = huxg_final;
results_constN.xgs = xgs;
results_constN.ts = ts;
results_constN.hs = hs';
results_constN.us = us';
results_constN.xg_f = xg_f;

%results_constN.time_to_ss = time_to_ss; 
fname = 'const_N_coulomb_5000.mat';
save(fname,'results_constN');

%% Implicit system of equations function (using discretization scheme from Schoof 2007)
function F = flowline_eqns(huxg,params)
    
    %vars unpack
    h = huxg(1:params.Nx);
    u = huxg(params.Nx+1:2*params.Nx);
    xg = huxg(2*params.Nx+1);
    hf = (-bed_schoof(xg.*params.x0,params)/params.h0)/(params.r);

    %grid params unpack
    dt = params.dt/params.t0;
    ds = params.dsigma;
    Nx = params.Nx;
    N1 = params.N1;
    sigma = params.sigma;
    sigma_elem = params.sigma_elem;
    b = -bed_schoof(xg.*sigma.*params.x0,params)/params.h0;
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
                  gamma.*N(1).* (u(1)./(u(1) + N(1).^nglen)).^m -...
                  0.5.*(h(1)+h(2)).*(h(2)-b(2)-h(1)+b(1))./(xg.*ds(1));
    Fu(2:Nx-1) = (alpha).*(1./(xg.*ds(2:Nx-1)).^((1/nglen)+1)).*...
                 (h(3:Nx).*(u(3:Nx)-u(2:Nx-1)).*abs(u(3:Nx)-u(2:Nx-1)).^((1/nglen)-1) -...
                  h(2:Nx-1).*(u(2:Nx-1)-u(1:Nx-2)).*abs(u(2:Nx-1)-u(1:Nx-2)).^((1/nglen)-1)) -...
                  gamma.*N(2:Nx-1).*(u(2:Nx-1)./(u(2:Nx-1) + N(2:Nx-1).^nglen)).^m -...
                  0.5.*(h(2:Nx-1)+h(3:Nx)).*(h(3:Nx)-b(3:Nx)-h(2:Nx-1)+b(2:Nx-1))./(xg.*ds(2:Nx-1));
    Fu(N1)     = (u(N1+1)-u(N1))/ds(N1) - (u(N1)-u(N1-1))/ds(N1-1);
    Fu(Nx)     = alpha.*(1./(xg.*ds(Nx-1)).^(1/nglen)).*...
                 (abs(u(Nx)-u(Nx-1)).^((1/nglen)-1)).*(u(Nx)-u(Nx-1)) - 0.5*(1-params.r)*hf;
             
    Fxg        = 3*h(Nx) - h(Nx-1) - 2*hf;   
    
    F = [Fh;Fu;Fxg];
end

