# livestock_disease_control_differing_social_perspectives

This repository contains files for performing computational simulations of a mathematical model framework to investigate the role of vaccine behaviour across a population of farmers on epidemic outbreaks amongst livestock, caused by pathogens with differential speed of spread over spatial landscapes of farms for two counties in England (Cumbria and Devon).

[![DOI](https://zenodo.org/badge/423575897.svg)](https://zenodo.org/badge/latestdoi/423575897)

The code was developed for the analysis presented in the scientific paper "Modelling livestock infectious disease control policy under differing social perspectives on vaccination behaviour" by Edward M. Hill, Naomi S. Prosser, Eamonn Ferguson, Jasmeet Kaler, Martin J. Green, Matt J. Keeling and Michael J. Tildesley.

Preprint details: EM Hill, NS Prosser, E Ferguson, J Kaler, MJ Green, MJ Keeling, MJ Tildesley. (2021) Modelling livestock infectious disease control policy under differing social perspectives on vaccination behaviour. agriRxiv. doi: 10.31220/agriRxiv.2021.00100. URL: https://doi.org/10.31220/agriRxiv.2021.00100

Model simulations are performed using the programming language Julia.
Julia makes use of environments, allowing bespoke package lists for separate projects. Documentation on working with environments and installing packages in the same state that is given by the project manifest: https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project

The main analyses in the paper require data that are not publicly available.  
We instead provide synthetic cattle and sheep demography data for Cumbria and Devon, with premises (matching the number retained in the data used for our analyses) randomly located in a 100kmx100km region. This enables the user to run the main simulation script (run_generic_livestock_disease_spatial_model_sellke_vers.jl) and a supporting script to compute the probability of a vaccinated premises not having been infected if it had been unvaccinated (ComputeVaccButNoInfecStats.jl). Note the results produced using the synthetic data would not match those from the paper.

We provide results files that allow the user to run the script returning the optimal strategy from a given perspective (ComputeOptimalRiskThreshold_VaccBehav_Script.jl) and the plotting scripts.

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
Main plotting scripts are produce_manuscript_plots.jl and produce_manuscript_plots.m.

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
