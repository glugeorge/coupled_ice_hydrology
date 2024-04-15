load('S1C_sensitivity_highres.mat');
%% Plotting

fig1 = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
M_vals = ["1.00" "4.64" "21.5" "100"];
A_vals = ["3.9" "8.6" "19" "42"];
t = tiledlayout(length(a_vals),length(C_vals),'TileSpacing','Compact');
for i=1:length(a_vals)
    for j=1:length(C_vals)
        nexttile;
        h = heatmap(A_vals,M_vals,squeeze(N_pos_arr(i,j,:,:))','ColorbarVisible','off','ColorLimits',[0.7 1],'Colormap',parula);
        %set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %ytickformat('%.2f');
        %xtickformat('%.2f');
         h.CellLabelFormat = '%.2f';
         if i == 4
             h.XLabel = ['C=',num2str(C_vals(j),'%.2f')];
         end
         if j == 1
             
             h.YLabel = ['a=',num2str(a_vals(i),'%.2f')];
         end
         %h.XDisplayLabels = {'SM','MED','LG'};
        %title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])
    end

end
ax = axes(t,'visible','off','Colormap',h.Colormap,'CLim',[0.7 1]);
title(t,'Sensitivity of peak in N');
cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.Label.String = 'relative position of maximum N';
%%
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
M_vals = ["1.00" "4.64" "21.5" "100"];
A_vals = ["3.9" "8.6" "19" "42"];
fig2 = figure(2);
t = tiledlayout(length(a_vals),length(C_vals),'TileSpacing','Compact');
for i=1:length(a_vals)
    for j=1:length(C_vals)
        nexttile;
        h = heatmap(A_vals,M_vals,squeeze(1e3.*xg_arr(i,j,:,:))','ColorbarVisible','off','ColorLimits',[30 500],'Colormap',parula);
        %set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %ytickformat('%.2f');
        %xtickformat('%.2f');
         h.CellLabelFormat = '%.0f';
         if i == 4
             h.XLabel = ['C=',num2str(C_vals(j),'%.2f')];
         end
         if j == 1
             
             h.YLabel = ['a=',num2str(a_vals(i),'%.2f')];
         end
         %h.XDisplayLabels = {'SM','MED','LG'};
        %title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])
    end

end
ax = axes(t,'visible','off','Colormap',h.Colormap,'CLim',[30 500]);
title(t,'Sensitivity of grounding line position')

cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.Label.String = 'grounding line position [km]';

%%
load('S1B_sensitivity.mat');

fig1 = figure(3);
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
M_vals = ["1.00" "4.64" "21.5" "100"];
A_vals = ["3.9" "8.6" "19" "42"];
t = tiledlayout(length(a_vals),1,'TileSpacing','Compact');
for i=1:length(a_vals)
        nexttile;
        h = heatmap(A_vals,M_vals,squeeze(N_pos_arr(i,:,:))','ColorbarVisible','off','ColorLimits',[0.4 1],'Colormap',parula);
        %set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %ytickformat('%.2f');
        %xtickformat('%.2f');
         h.CellLabelFormat = '%.2f';
             
         h.YLabel = ['a=',num2str(a_vals(i),'%.2f')];
         if i ==4
             h.XLabel = ''
         end
         %h.XDisplayLabels = {'SM','MED','LG'};
        %title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])

end
ax = axes(t,'visible','off','Colormap',h.Colormap,'CLim',[0.4 1]);
title(t,'Sensitivity of peak in N');
cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.Label.String = 'relative position of maximum N';
%%
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
M_vals = ["1.00" "4.64" "21.5" "100"];
A_vals = ["3.9" "8.6" "19" "42"];
fig2 = figure(4);
t = tiledlayout(length(a_vals),1,'TileSpacing','Compact');
for i=1:length(a_vals)
        nexttile;
        h = heatmap(A_vals,M_vals,squeeze(1e3.*xg_arr(i,:,:))','ColorbarVisible','off','ColorLimits',[30 350],'Colormap',parula);
        %set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %ytickformat('%.2f');
        %xtickformat('%.2f');
         h.CellLabelFormat = '%.0f';
             
         h.YLabel = ['a=',num2str(a_vals(i),'%.2f')];
         %h.XDisplayLabels = {'SM','MED','LG'};
        %title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])

end
ax = axes(t,'visible','off','Colormap',h.Colormap,'CLim',[30 350]);
title(t,'Sensitivity of grounding line position')

cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.Label.String = 'grounding line position [km]';
