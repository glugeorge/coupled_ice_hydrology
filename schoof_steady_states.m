% Schoof steady states

%% Physical parameters
params.n = 3;
params.rho_i = 917;
params.rho_w = 1028;
params.g = 9.81;
params.C = 0.5; % base 0.2
params.As = 2.26e-21; % Calculated 
params.f = 0.07; % From Kingslake thesis
params.K0 = 10^-24; % From Kingslake thesis  
params.L = 3.3e5; % Kingslake thesis
params.year = 3600*24*365;
%% Scaling params (coupled model equations solved in non-dim form)
params.x0 = 1000*10^3;
params.h0 = 1000;
params.Q0 = 1500;

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
%gamma = As*(C*N0)^n;
params.gamma = (params.C*params.N0*params.x0)/(params.rho_i*params.g*params.h0^2);
params.beta = params.th0/params.t0;
params.r = params.rho_i/params.rho_w;

params.transient = 0;

%% Grid parameters - ice sheet
params.Nx = 600;                    %number of grid points - 200
params.N1 = 100;                    %number of grid points in coarse domain - 100
params.Nh = 600;
params.sigGZ = 0.85;                %extent of coarse grid (where GL is at sigma=1) - 0.97
sigma1=linspace(params.sigGZ/(params.N1+0.5), params.sigGZ, params.N1);
sigma2=linspace(params.sigGZ, 1, params.Nx-params.N1+1);
params.sigma = [sigma1, sigma2(2:end)]';    %grid points on velocity (includes GL, not ice divide)
params.dsigma = diff(params.sigma); %grid spacing
params.sigma_elem = [0;(params.sigma(1:params.Nx-1) + params.sigma(2:params.Nx))./2 ]; %grid points on thickness (includes divide, not GL)

%% Grid parameters - hydro
params.sigma_h = linspace(0,1,params.Nh)';
params.dsigma_h = diff(params.sigma_h); %grid spacing

%% Determine at what points there is coupling
% 1 - coupling on, 0 - coupling off
params.hydro_u_from_ice_u = 1;
params.hydro_psi_from_ice_h = 1;
params.ice_N_from_hydro = 1;

%% Initial Guess conditions
params.shear_scale = 1;
Q = 0.001*ones(params.Nh,1);
N = ones(params.Nh,1);
S = 5/params.S0*ones(params.Nh,1); 
params.S_old = S;
params.M = 0e-4/params.M0; % zero when using schoof bed
params.N_terminus = 0;
params.accum = 1./params.year;
xg = 1500e3/params.x0; % Set high past sill for retreat
hf = (-bed_schoof(xg.*params.x0,params)/params.h0)/params.r;
h = 1 - (1-hf).*params.sigma;
u = 0.1*(params.sigma_elem.^(1/3)) + 1e-3; % 0.1 for C = 0.5, 0.3 for C = 0.1-0.4
params.Q_in = 10/params.Q0;

params.h_old = h;
params.xg_old = xg;
params.ice_start = 3*params.Nh;

sig_old = params.sigma;
sige_old = params.sigma_elem;
QNShuxg0 = [Q;N;S;h;u;xg];
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 
params.dt = 1; % arbitrary since steady state

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
for A=[0.9e-25, 1.9e-25]%, 2.9e-25, 3.9e-25, 4.9e-25]
    params.A = A;
    params.alpha = 2*params.u0^(1/params.n)/(params.rho_i*params.g*params.h0*(params.x0*params.A)^(1/params.n));

    flf = @(QNShuxg) schoof_combined_hydro_ice_eqns(QNShuxg,params);
    
    [QNShuxg_init,F,exitflag,output,JAC] = fsolve(flf,QNShuxg0,options);
    
    Q = QNShuxg_init(1:params.Nh);
    N = QNShuxg_init(params.Nh+1:2*params.Nh);
    S = QNShuxg_init(2*params.Nh+1:3*params.Nh);
    h = QNShuxg_init(params.ice_start+1:params.ice_start+ params.Nx);
    u = QNShuxg_init(params.ice_start + params.Nx+1:params.ice_start+2*params.Nx);
    xg = QNShuxg_init(params.ice_start+2*params.Nx+1);

    ax1 = subplot(5,1,1); 
    plot(sigma_h.*xg.*params.x0./1000,Q); ylabel('Q');title('Nondimensionalized steady state variables');
    hold on;
    text(max(sigma_h.*xg.*params.x0./1000),max(Q),['A=',num2str(A*1e25)])

    ax2 =subplot(5,1,2);
    plot(sigma_h.*xg.*params.x0./1000,N);ylabel('N');hold on;
    ax3 = subplot(5,1,3);
    plot(sigma_h.*xg.*params.x0./1000,S);ylabel('S');hold on;
    ax4 = subplot(5,1,4);
    plot(sigma_elem.*xg.*params.x0./1000,h);ylabel('h');hold on;
    ax5 = subplot(5,1,5);
    plot(sigma.*xg.*params.x0./1000,u);hold on;ylabel('u');xlabel('distance from divide, \emph{x} (km)','Interpreter','latex')
    linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
    xlim([0,1500]);

end
subplot(5,1,4)
plot(linspace(0,1)*1500,bed_schoof(linspace(0,1).*1500e3,params)./params.h0,'-k'); 
