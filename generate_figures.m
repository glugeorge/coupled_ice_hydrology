%% This script generates the figures for Lu & Kingslake, in prep
% Running will create all figures as once

%% Figure 1
% This figure is a schematic drawn in powerpoint. 

%% Figure 2
run 'coulomb_sensitivity/plot_sensitivity.m';

%% Figure 3
run 'budd_sensitivity/plot_sensitivity.m';

%% Figure 4
run 'results_smallQin/steady_state.m';

%% Figure 5
run 'results_smallQin/term_plot.m';

%% Figure 6
run 'results_smallQin/trans_evolution.m';

%% Figure 7
run 'results_smallQin/compare_transient.m';
