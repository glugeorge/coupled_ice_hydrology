bed_func = @bed_schoof;

%% Physical parameters
params.Cf = 0.4;
params.A = 1e-25; 
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
params.h0 = 1000;
params.Q0 = 1;

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
params.Nx = 700;                    %number of grid points - 200
params.N1 = 100;                    %number of grid points in coarse domain - 100
params.sigGZ = 0.85;                %extent of coarse grid (where GL is at sigma=1) - 0.97
sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
params.sigma = [sigma1, sigma2(2:end)]';    %grid points on velocity (includes GL, not ice divide)
params.dsigma = diff(params.sigma); %grid spacing
params.sigma_elem = [0;(params.sigma(1:params.Nx-1) + params.sigma(2:params.Nx))./2 ]; %grid points on thickness (includes divide, not GL)


%% Establish timings
params.year = 3600*24*365;  %number of seconds in a year
params.Nt =1000;                    %number of time steps - normally 150
params.end_year = 1000; %normally 7500

params.dt = params.end_year*params.year/params.Nt;

%% Determine at what points there is coupling
% 1 - coupling on, 0 - coupling off
params.hydro_u_from_ice_u = 1;
params.hydro_psi_from_ice_h = 1;
params.ice_N_from_hydro = 1;

%% Initial "steady state" conditions
params.shear_scale = 1;
params.accum = 1./params.year;
xg = 1500e3/params.x0; % Set high past sill for retreat
hf = (-bed_func(xg.*params.x0,params)/params.h0)/params.r;
h = 1 - (1-hf).*params.sigma;
u = 0.1*(params.sigma_elem.^(1/3)) + 1e-3; % 0.1 for C = 0.5, 0.3 for C = 0.1-0.4
params.Q_in = 0.01/params.Q0;

params.h_old = h;
params.xg_old = xg;

sig_old = params.sigma;
sige_old = params.sigma_elem;
huxg0 = [h;u;xg];

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
flf = @(huxg) flowline_eqns(huxg,params,bed_func);

[huxg_init,F,exitflag,output,JAC] = fsolve(flf,huxg0,options);


h = huxg_init(1:params.Nx);
u = huxg_init(params.Nx+1:2*params.Nx);
xg = huxg_init(2*params.Nx+1);
hf = (-bed_func(xg.*params.x0,params)/params.h0)/(params.r);


%% Calculate transient GL evolution over bedrock peak
huxg_t = huxg_init;
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
C_fs = linspace(0.4,1,10);
for t=1:params.Nt
    if t <= 10
        params.Cf = C_fs(t);
    else
        params.Cf = C_fs(10);
    end
    flf = @(huxg) flowline_eqns(huxg,params,bed_func);
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
load coulomb_retreat.mat
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
results_simpleN.params = params;
results_simpleN.init_cond = huxg_init;
%results_constN.steady_state = huxg_final;
results_simpleN.xgs = xgs;
results_simpleN.ts = ts;
results_simpleN.hs = hs';
results_simpleN.us = us';
%results_constN.xg_f = xg_f;

%results_constN.time_to_ss = time_to_ss; 
fname = 'simple_N_coulomb.mat';
save(fname,'results_simpleN');

%% Implicit system of equations function (using discretization scheme from Schoof 2007)
function F = flowline_eqns(huxg,params,bed_func)
    
    %vars unpack
    h = huxg(1:params.Nx);
    u = huxg(params.Nx+1:2*params.Nx);
    xg = huxg(2*params.Nx+1);
    hf = (-bed_func(xg.*params.x0,params)/params.h0)/(params.r);

    %grid params unpack
    dt = params.dt/params.t0;
    ds = params.dsigma;
    Nx = params.Nx;
    N1 = params.N1;
    sigma = params.sigma;
    sigma_elem = params.sigma_elem;
    b = -bed_func(xg.*sigma.*params.x0,params)/params.h0;
    Fh = zeros(Nx,1);
    Fu = zeros(Nx,1);
    Fn = zeros(Nx,1);

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
    
    h_interp = interp1(sigma_elem,h.*params.h0,sigma,'linear','extrap');
    N  = (params.g.*params.rho_i.*h_interp(1:end) - max(params.rho_w.*params.g.*params.h0.*b(1:end),zeros(Nx,1)))./params.N0;
    

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
                 (abs(u(Nx)-u(Nx-1)).^((1/nglen)-1)).*(u(Nx)-u(Nx-1)) - params.Cf.*0.5*(1-params.r)*hf;
             
    Fxg        = 3*h(Nx) - h(Nx-1) - 2*hf;  

   
    F = [Fh;Fu;Fxg];
end

