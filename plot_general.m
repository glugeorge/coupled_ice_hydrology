folder = 'C_sensitivity';
myFiles = dir(fullfile(folder,'*.mat'));

for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(folder, baseFileName);
  load(fullFileName);
  if mod(k,2) == 0
      hold on;
      % grounding line evolution
      plot(results_constN.ts,results_constN.xgs.*results_constN.params.x0./1e3,'linewidth',3,'DisplayName','Constant N');xlabel('time (yr)');ylabel('x_g');
      legend;
  else
      subplot(length(myFiles)/4,2,floor(k/2)+1);
      % grounding line evolution
      plot(results.ts,results.xgs.*results.params.x0./1e3,'linewidth',3,'DisplayName','Coupled');xlabel('time (yr)');ylabel('x_g');
      new_title = replace(baseFileName,'_A',', A');
      new_title = replace(new_title,'_c.mat','');
      new_title = replace(new_title,'_','=');
      title(new_title);
  end
end
