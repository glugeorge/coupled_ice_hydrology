year = 60*60*24*365;
n = 3;
figure();
u_arr = linspace(1/year,2000/year);
u_mesh = ones(length(u_arr)).*u_arr;
N_arr = linspace(0.01,300000);
N_mesh = (ones(length(N_arr)).*N_arr)';
C = [0.01 0.1 1];
As = [2.26e-22,2.26e-21,2.26e-20];
c_limit = [0,3e5];
for i = 1:3
    for j = 1:3
        subplot(3,3,(i-1)*3 + j);
        tau = C(i)*N_mesh.*(u_mesh./(u_mesh + As(j)*(C(i)*N_mesh).^n)).^(1/n);
        surface(u_mesh*year,N_mesh,tau,'EdgeColor','none');colorbar;xlabel('u (m/yr)');ylabel('N (Pa)');
        %caxis(c_limit);set(gca,'ColorScale','log')
        title(['C=',num2str(C(i)),'; As=',num2str(As(j))]);

    end
end
sgtitle('Basal shear stress as a function of N, u, As, C');

