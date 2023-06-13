%% Plotting
close all; clear;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
myFiles = dir(fullfile('*.mat'));
% first load variables to plot
figure(Position=[250,550,800,200])
tiledlayout('horizontal',"Padding","loose",'TileSpacing','tight')
for i=1:length(myFiles)
    load(myFiles(i).name);
end
variables = [A_arr*1e25; accums; M_ins*1e4; Q_ins; As_ar*1e21; C_arr];
xlabels = ["$A$ [$\times10^{-25}$ s$^{-1}$ Pa$^{-3}$]",...
            "$a$ [m yr$^{-1}$]",...
            "$M$ [$\times10^{-4}$ m$^2$ s$^{-1}$]",...
            "$Q_{\mathrm{in}}$ [m$^3$ s$^{-1}$]",...
            "$A_s$ [$\times10^{-21}$ m s$^{-1}$ Pa$^{-3}$]",...
            "$C_C$"]; 
titles = [""];
values_used = [2.9,1,0,10,2.26,0.2]
for i=1:length(myFiles)
    load(myFiles(i).name);
    to_plot = variables(i,:);
    nexttile;
    contourf(to_plot,params.sigma_elem,hs'.*params.h0/1000,EdgeColor='None'); 
    caxis([0 2]);
    xlabel(xlabels(i),'Interpreter','latex');
    set(gca,'Ydir','Reverse');
    xline(values_used(i), '--k',LineWidth=2);
    ylim([0.5 1])
    if i == 1
        ylabel('normalized distance from divide, $\sigma$','Interpreter','latex')
    else
        set(gca,'YTickLabel',[]);
    end
end
c = colorbar;
c.Label.String = 'ice thickness $h$ [km]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

figure(Position=[250,550,800,200])
tiledlayout('horizontal',"Padding","loose",'TileSpacing','tight')
for i=1:length(myFiles)
    load(myFiles(i).name);
    to_plot = variables(i,:);
    nexttile;
    contourf(to_plot,params.sigma_h,Ns'.*params.N0./1e6,EdgeColor='None'); 
    caxis([0 1]);
    xlabel(xlabels(i),'Interpreter','latex');
    xline(values_used(i), '--k',LineWidth=2);
    set(gca,'Ydir','Reverse');
    ylim([0.5 1])
    if i == 1
        ylabel('normalized distance from divide, $\sigma$','Interpreter','latex')
    else
        set(gca,'YTickLabel',[]);
    end
end
c = colorbar;
c.Label.String = 'effective pressure $N$ [MPa]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

% tiledlayout('horizontal',"Padding","loose",'TileSpacing','tight');
% load sens_C.mat
% to_plot = C_arr;%As_ar*1e21; %Q_ins%M_ins*1e4;%accums;%A_arr*1e25;
%x_label = '$A$ [$\times10^{-25}$ s$^{-1}$ Pa$^{-3}$]';
%x_label = '$a$ [m yr$^{-1}$]';
%x_label = '$M$ [$\times10^{-4}$ m$^2$ s$^{-1}$]';
%x_label = '$Q_{\mathrm{in}}$ [m$^3$ s$^{-1}$]';
%x_label = '$A_s$ [$\times10^{-21}$ m s$^{-1}$ Pa$^{-3}$]';
% x_label = '$C_C$';
% 
% a1 = nexttile;
% contourf(to_plot,params.sigma_elem,hs'.*params.h0); hold on;
% caxis([0 2500]);
% 
% set(gca,'Ydir','Reverse');
% ylim([0.5 1]);
% xlabel(x_label,Interpreter='latex');
% ylabel('Normalized distance from divide',Interpreter='latex');
% title('Ice Thickness',Interpreter='latex')

% a2 = nexttile;
% contourf(to_plot,params.sigma_h,Ns'.*params.N0./1e6,EdgeColor='None',FaceAlpha=0.75);hold on;
% [C,h]=contour(to_plot,params.sigma_h,Ns'.*params.N0./1e6,"LabelFormat","%0.1f MPa",'LabelSpacing',200);
% %clabel(C,h,[0.5 0.6 0.7 0.8 0.9 1]);
% caxis([0 1]);
% 
% ylim([0.5 1]);
% xlabel(x_label,Interpreter='latex');
% ylabel('Normalized distance from divide',Interpreter='latex');
% title('Effective Pressure',Interpreter='latex')
% set(gca,'Ydir','Reverse');
% linkaxes([a1,a2],'x');

%% accum
% load sens_accum.mat
% a1 = nexttile;
% contourf(accums,params.sigma_elem,hs'.*params.h0,EdgeColor='None',FaceAlpha=0.75); hold on;
% [C,h]=contour(A_arr,params.sigma_elem,hs'.*params.h0,"LabelFormat","%d m");
% clabel(C,h);
% 
% set(gca,'Ydir','Reverse');
% set(gca,'XTickLabel','')
% ylim([0.5 1]);
% a2 = nexttile;
% contourf(accums,params.sigma_h,Ns'.*params.N0./1e6,EdgeColor='None',FaceAlpha=0.75);hold on;
% [C,h]=contour(A_arr,params.sigma_h,Ns'.*params.N0./1e6,"LabelFormat","%0.1f MPa");
% clabel(C,h,[0.5 0.8 1]);
% 
% ylim([0.5 1]);
% %set(gca,'Ydir','Reverse');
% linkaxes([a1,a2],'x');
