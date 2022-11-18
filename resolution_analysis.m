folder = 'results/spatial_coords/hydro';
myFiles = dir(fullfile(folder,'*.mat'));
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(folder, baseFileName);
  load(fullFileName);
  Nx = results.params.Nx;
  N1 = results.params.N1;
  Nh = results.params.Nh;
  linename = strcat('Nh=',num2str(Nh));
  split = results.params.sigGZ;
  subplot(1,2,1);
  plot(results.ts.*results.params.t0,results.xgs.*results.params.x0./1000,'o','DisplayName',linename);
  hold on;
  subplot(1,2,2);
  QNShuxg = results.steady_state;
  %h = QNShuxg(3*Nx+1:4*Nx);
  %xg = QNShuxg(5*Nx+1);
  %N = QNShuxg(Nx+1:2*Nx);
  h = QNShuxg(3*Nh+1:3*Nh+1+Nx);
  xg = QNShuxg(3*Nh+2*Nx+1);
  N = QNShuxg(Nh+1:2*Nh);
  plot(results.params.sigma_h.*xg.*results.params.x0./1000,N,'DisplayName',linename)
  hold on;
end
subplot(1,2,1);
legend;
ylabel('grounding line position, \emph{$x_g$} (km)','Interpreter','latex');
xlabel('time, \emph{t} (y)','Interpreter','latex');
title('Grounding line evolution');

subplot(1,2,2);
legend;
ylabel('ice thickness, \emph{h} (m)','Interpreter','latex');
xlabel('distance from divide, \emph{x} (km)','Interpreter','latex');
title('Final steady state profile'); 
