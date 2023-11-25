# coupled_ice_hydrology
## A code repository for: "Coupling between ice flow and subglacial hydrology enhances marine ice-sheet retreat"
## George Lu and Jonathan Kingslake

This repository holds the run scripts for the model described by Lu and Kingslake (in prep). 

To generate the foundational figures for Lu and Kingslake (in prep), open MATLAB, and run
`generate_figures.m`.

To run the experiments described in Lu and Kingslake (in prep), follow the following scripts:
Sensitivity: To test the sensibility to a parameter, go to the either the `budd_sensitivity` or `coulomb_sensitivity` folder, and select the the script `sensitivity_PARAM.m`. 

Experiment S1.B: `steady_budd.m`
Experiment S1.C: `steady_coulomb.m`
Experiment S2: `oneway_hydro.m`

Experiment T1.B: `transient_budd.m`
Experiment T1.C: `transient_coulomb.m`
Experiment T2.B: `staticN_transient_budd.m`
Experiment T2.C: `staticN_transient_coulomb.m`