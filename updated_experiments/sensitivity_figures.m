load('S1C_sensitivity.mat');
%% Plotting

fig1 = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
M_vals = ["1.00" "4.64" "21.5" "100"];
A_vals = ["3.9" "8.6" "19" "42"];
t = tiledlayout(length(a_vals),length(C_vals),'TileSpacing','Compact');
for i=1:length(a_vals)
    for j=1:length(C_vals)
        nexttile;
        arr = squeeze(N_pos_arr(i,j,:,:))';
        if i == 4 && j==1
            h = heatmap(A_vals,M_vals,arr,'ColorbarVisible','off','ColorLimits',[0.5 1],'Colormap',parula);
        else
            h = heatmap(arr,'XDisplayLabels',NaN*ones(length(arr),1),'YDisplayLabels',NaN*ones(length(arr),1),'ColorbarVisible','off','ColorLimits',[0.5 1],'Colormap',parula);
        end%set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %ytickformat('%.2f');
        %xtickformat('%.2f');
         h.CellLabelFormat = '%.2f';
         h.Interpreter= 'latex';
         %if i == 4
         %    h.XLabel = ['C=',num2str(C_vals(j),'%.2f')];
             
         %end
         %if j == 1
             
         %    h.YLabel = ['a=',num2str(a_vals(i),'%.2f')];
         %end
         %h.XDisplayLabels = {'SM','MED','LG'};
        %title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])
    end

end
ax = axes(t,'visible','off','Colormap',h.Colormap,'CLim',[0.5 1]);
title(t,'Sensitivity of peak in N','Interpreter','latex');
cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.Label.String = 'relative position of maximum N';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
%%
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
M_vals = ["1.00" "4.64" "21.5" "100"];
A_vals = ["3.9" "8.6" "19" "42"];
fig2 = figure(2);
t = tiledlayout(length(a_vals),length(C_vals),'TileSpacing','Compact');
for i=1:length(a_vals)
    for j=1:length(C_vals)
        nexttile;
        arr = squeeze(1e2.*xg_arr(i,j,:,:))';
        if i == 4 && j==1
            h = heatmap(A_vals,M_vals,arr,'ColorbarVisible','off','ColorLimits',[30 500],'Colormap',parula);
        else
            h = heatmap(arr,'XDisplayLabels',NaN*ones(length(arr),1),'YDisplayLabels',NaN*ones(length(arr),1),'ColorbarVisible','off','ColorLimits',[30 500],'Colormap',parula);
        end
        %set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %ytickformat('%.2f');
        %xtickformat('%.2f');
         h.CellLabelFormat = '%.0f';
         h.Interpreter= 'latex';
         % if i == 4
         %     h.XLabel = ['C=',num2str(C_vals(j),'%.2f')];
         % end
         % if j == 1
         % 
         %     h.YLabel = ['a=',num2str(a_vals(i),'%.2f')];
         % end
         %h.XDisplayLabels = {'SM','MED','LG'};
        %title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])
    end

end
ax = axes(t,'visible','off','Colormap',h.Colormap,'CLim',[30 500]);
title(t,'Sensitivity of grounding line position','Interpreter','latex')

cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.Label.String = 'grounding line position [km]';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';

%%
load('S1B_sensitivity.mat');

fig1 = figure(3);
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
M_vals = ["1.00" "4.64" "21.5" "100"];
A_vals = ["3.9" "8.6" "19" "42"];
t = tiledlayout(length(a_vals),1,'TileSpacing','Compact');
for i=1:length(a_vals)
        nexttile;
        arr = squeeze(N_pos_arr(i,:,:))';
        if i == 4
            h = heatmap(A_vals,M_vals,arr,'ColorbarVisible','off','ColorLimits',[0.5 1],'Colormap',parula);
        else
            h = heatmap(arr,'XDisplayLabels',NaN*ones(length(arr),1),'YDisplayLabels',NaN*ones(length(arr),1),'ColorbarVisible','off','ColorLimits',[0.5 1],'Colormap',parula);
        end%set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %ytickformat('%.2f');
        %xtickformat('%.2f');
         h.CellLabelFormat = '%.2f';
         h.Interpreter= 'latex';
         %h.YLabel = ['a=',num2str(a_vals(i),'%.2f')];
         if i ==4
             h.XLabel = ''
         end
         %h.XDisplayLabels = {'SM','MED','LG'};
        %title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])

end
ax = axes(t,'visible','off','Colormap',h.Colormap,'CLim',[0.5 1]);
title(t,'Sensitivity of peak in N','Interpreter','latex');
cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.Label.String = 'relative position of maximum N';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
%%
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
M_vals = ["1.00" "4.64" "21.5" "100"];
A_vals = ["3.9" "8.6" "19" "42"];
fig2 = figure(4);
t = tiledlayout(length(a_vals),1,'TileSpacing','Compact');
for i=1:length(a_vals)
        nexttile;
        arr = squeeze(1e2.*xg_arr(i,:,:))';
        if i == 4
            h = heatmap(A_vals,M_vals,arr,'ColorbarVisible','off','ColorLimits',[30 500],'Colormap',parula);
        else
            h = heatmap(arr,'XDisplayLabels',NaN*ones(length(arr),1),'YDisplayLabels',NaN*ones(length(arr),1),'ColorbarVisible','off','ColorLimits',[30 500],'Colormap',parula);
        end%set(gca, 'XScale', 'log');
        %set(gca, 'YScale', 'log');
        %ytickformat('%.2f');
        %xtickformat('%.2f');
         h.CellLabelFormat = '%.0f';
         h.Interpreter= 'latex';    
         h.YLabel = ['a=',num2str(a_vals(i),'%.2f')];
         %h.XDisplayLabels = {'SM','MED','LG'};
        %title(['a=',num2str(a_vals(i)),', C=',num2str(C_vals(j))])

end
ax = axes(t,'visible','off','Colormap',h.Colormap,'CLim',[30 500]);
title(t,'Sensitivity of grounding line position','Interpreter','latex')

cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.Label.String = 'grounding line position [km]';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';