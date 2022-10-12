%% Case comparing with power law
load coupled_init_cond.mat;
n = params.n;
m = 1/n;
C2 = params.C;
N = interp1(params.sigma_h,params.N0.*N,params.sigma,"linear","extrap");
u = u.*params.u0;
new_tau = C2.*N.*(u./(u+params.As*(C2*N).^n)).^m;
C_new = new_tau./(abs(u).^(m-1).*u);

%% Case comparing with constant N coulomb friction law
params.const_N = ones(size(N)).*100000;
params.new_tau = new_tau;
params.u = u;
C_init = C2.*ones(size(params.const_N));
flf = @(C) solve_nonuni_C(C,params);
options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);

[C_new,F,exitflag,output,JAC] = fsolve(flf,C_init,options);

function F=solve_nonuni_C(C,params)
    F = params.new_tau - C.*params.const_N.*(params.u./(params.u+params.As.*(C.*params.const_N).^params.n)).^(1/params.n);

end