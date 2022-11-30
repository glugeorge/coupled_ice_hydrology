load C_0.5_A_2.9_c_schoof_retreat.mat;
params = results.params;
time = 1;
sigma_elem = params.sigma_elem;
sigma_h = params.sigma_h;
sigma = params.sigma; 
Q = results.Qs(:,time);
N = results.Ns(:,time);
S = results.Ss(:,time);
h = results.hs(:,time);
u = results.us(:,time);
xg = results.xgs(time);

N_ice = interp1(sigma_h,N,sigma);
max_N = max(N_ice);
mid_N = max_N/2;
max_u = max(u);
min_u = min(u);
figure();
plot(sigma,N_ice,'DisplayName','N');hold on; plot(sigma,u,'DisplayName','u');
title('Nondimensionalized variable steady state profiles');
xlabel('sigma'); legend;

figure();
u_arr = linspace(min_u,max_u);
tau = params.gamma.*0.01.*(u_arr./(u_arr + 0.01.^params.n)).^(1/params.n);
plot(u_arr,tau,'DisplayName',['N = ',num2str(0.01)]);
hold on;
tau = params.gamma.*mid_N.*(u_arr./(u_arr + mid_N.^params.n)).^(1/params.n);
plot(u_arr,tau,'DisplayName',['N = ',num2str(mid_N)]);

tau = params.gamma.*max_N.*(u_arr./(u_arr + max_N.^params.n)).^(1/params.n);
plot(u_arr,tau,'DisplayName',['N = ',num2str(max_N)]);
ylabel('Non-dimensionalized tau');xlabel('Non-dimensionalized u');legend