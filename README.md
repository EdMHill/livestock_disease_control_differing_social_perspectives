# livestock_disease_control_differing_social_perspectives

This repository contains files for performing computational simulations of a mathematical model framework to investigate the role of vaccine behaviour across a population of farmers on epidemic outbreaks amongst livestock, caused by pathogens with differential speed of spread over spatial landscapes of farms for two counties in England (Cumbria and Devon).

The code was developed for the analysis presented in the scientific paper "Modelling disease control policy under differing social perspectives: An application to infectious diseases of livestock" by Edward M. Hill, Naomi S. Prosser, Eamonn Ferguson, Jasmeet Kaler, Martin J. Green, Matt J. Keeling and Michael J. Tildesley.

Preprint details: TBC

Model simulations are performed using the programming language Julia.
Julia makes use of environments, allowing bespoke package lists for separate projects. Documentation on working with environments and installing packages in the same state that is given by the project manifest: https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project

The main analyses in the paper require data not publicly available. To enable the main simulations scripts to be run, we provide synthetic livestock demography and premises location data.

Please find below an explainer of the directory structure within this repository.

## data

**dummy_data**  
Directory containing synthetic data on livestock counts and premises location data.

**Shapefiles**  
County and country boundary shapefiles for UK.

## results
Directories to store simulation outputs and plot scripts.

**GB_county_model_simn_outputs**  
Simulation output files for Cumbria and Devon configurations.

**generate_plot_scripts**  
Plotting scripts

**GenericLandscapeSimnOutputs**  
Save location for output files from runs using dummy data.

## src

**matlab_packages**  
Supporting packages used to produce the ternary plots in MATLAB.

**spatial_simns**  
Houses directories with optimal strategy calculation code and mathematical model simulation files.

- **optim_behaviour_generic_model**
    * ComputeVaccButNoInfecStats.jl: Compute probability of vaccinated premises not being infected (after moment of vaccination).
    * ComputeOptimalRiskThreshold_VaccBehav_Script.jl: Find optimal risk threshold to minimise cost.

- **run_model_scripts**
    * run_generic_livestock_disease_spatial_model_sellke_vers.jl: Run the individual-based livestock disease epidemic model using the Sellke construction.

- **sellke_simn_fns**
    * Collection of files that contain the functions used in the main outbreak simulation.
