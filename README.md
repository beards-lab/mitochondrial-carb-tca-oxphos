# mitochondrial-carb-tca-oxphos
Computational model of mitochondrial metabolism integrating TCA cycle, carbohydrate transport, and oxidative phosphorylation.

This package simulates the mitochondrial metabolism model described in "Systems Analysis of Carboxylate Transport and Oxidation Pathways in Cardiac Mitochondria" by Collins et al. The codes are in MATLAB. Data are stored in Excel spreadsheets.

dXdT.m: This script computes the rates of change of the state variables in the model (the right-hand side of the ODE system). 

maincode.m: This script initiates the model and runs and plots the outputs of the kinetic time-course experiments.

Load_Data.m: This script loads the data and the substrate concentrations for the simulation experiments.

dXdT_electrode.m: This is the ODE model for the first-order filter used to simulate the electrode response for the oxygen consumption measurements.

Quasi_Steady_Experiments.m: This script simulates the quasi-steady experiments.

Substrates.xlsx: This spreadsheet includes all of the substrate concentrations used for the simulated experiments.

Metabolite_Data.xlsx: This spreadsheet includes the measured metabolic concentrations.

NADH_Data.xlsx: This spreadsheet includes all of the NAD(P)H measurements for the kinetic experiments.

Oxygen_Data.xlsx: This spreadsheet includes all of the oxygen data for the kinetic experiments.
