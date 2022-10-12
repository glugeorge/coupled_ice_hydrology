load coupled_init_cond.mat
plot(params.sigma_elem*xg*params.x0/1000,h*params.h0,'DisplayName','Coupled');
hold on;

load oneway_init_cond.mat
plot(params.sigma_elem*xg*params.x0/1000,h*params.h0,'o','DisplayName','Uncoupled');
