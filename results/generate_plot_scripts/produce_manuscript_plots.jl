#=
Purpose:
Script to call functions producing figures for the manuscript YYY

Useful Julia plotting resource:
https://nextjournal.com/leandromartinez98/tips-to-create-beautiful-publication-quality-plots-in-julia

Font family information:
https://gr-framework.org/fonts.html

Julia version: 1.6.3
Date: 3rd November 2021
=#

#--------------------------------------------------------------------------
# LOAD ENVIRONMENT
#--------------------------------------------------------------------------

#Set path to directory this file resides in
cd(dirname(@__FILE__))

using Pkg
Pkg.activate("../../../")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
using DelimitedFiles
using Plots, StatsPlots
using Plots.PlotMeasures # to use measures explicitly, e.g. mm
using TernaryPlots
using LaTeXStrings
using ColorSchemes
using MAT
using StatsBase
using Shapefile

#-------------------------------------------------------------------------------
# IMPORT REQUIRED FUNCTION FILES
#-------------------------------------------------------------------------------
include("OptimThresholdVsCostRatio_SellkeRuns_GenerateLinePlots.jl")
include("../../../Data/DataVis/2020_livestock_popn_data_plots/county_comparison/between_county_livestock_comparison_script.jl")
include("OptimThresholdVsCostRatio_SellkeRuns_GenerateTernaryPlots.jl")

#-------------------------------------------------------------------------------
# DEFINE RECIPIES FOR USE IN THE PLOTS
#-------------------------------------------------------------------------------
# Set spacing between markers
# Example usage: plot(range(0, 2π, length = 100), [sin, cos], st = :samplemarkers, step = 20, shape = :star)
@recipe function f(::Type{Val{:samplemarkers}}, x, y, z; step = 10)
    n = length(y)
    sx, sy = x[1:step:n], y[1:step:n]
    # add an empty series with the correct type for legend markers
    @series begin
        seriestype := :path
        markershape --> :auto
        x := []
        y := []
    end
    # add a series for the line
    @series begin
        primary := false # no legend entry
        markershape := :none # ensure no markers
        seriestype := :path
        seriescolor := get(plotattributes, :seriescolor, :auto)
        x := x
        y := y
    end
    # return  a series for the sampled markers
    primary := false
    seriestype := :scatter
    markershape --> :auto
    x := sx
    y := sy
end

#-------------------------------------------------------------------------------
# Cumbria & Devon livestock population demographics (Figure 1)
#-------------------------------------------------------------------------------

# Specify if generated figures should be saved to file
save_fig_flag = false

# Set county_IDs as a vector
county_IDs = [8,10]
# county_IDs = [0,0]
            # Cumbria: 8
            # Devon: 10
            # Dummy data: 0

# Call the function to produce the plot set
between_county_livestock_comparison_plots(save_fig_flag,county_IDs)

#-------------------------------------------------------------------------------
# All reactionary vaccinators scenario plots (Figure 2)
#-------------------------------------------------------------------------------

#Set path to directory this file resides in
cd(dirname(@__FILE__))

# Specify if generated figures should be saved to file
save_fig_flag = false

# Specify the configuration (location & pathogen attribute set) in use
# config_ID = "cumbria"
# config_ID = "cumbria_alt"
# config_ID = "devon"
# config_ID = "devon_alt"

# Call the function to produce the plot set
OptimThresholdVsCostRatio_SellkeRuns_GenerateLinePlots(save_fig_flag,config_ID)

#-------------------------------------------------------------------------------
# Optimal strategy ternary plots (Figure 3)
#-------------------------------------------------------------------------------

# Specify if generated figures should be saved to file
save_fig_flag = false

# Specify the configuration (location & pathogen attribute set) in use
config_ID = "cumbria"
# config_ID = "cumbria_alt"
# config_ID = "devon"
# config_ID = "devon_alt"

# Set scenario IDs & batch_ID_offset to be used
if (config_ID == "cumbria") || (config_ID == "devon")
    scen_IDs = collect(901:1:1131)
    batch_ID_offset = 15000
elseif (config_ID == "cumbria_alt") || (config_ID == "devon_alt")
    scen_IDs = collect(1201:1:1431)
    batch_ID_offset = 20000
else
    error("Invalid config_ID provided.")
end

# Relative cost of vaccination to generate plots for
rel_cost_of_vacc_idxs = [2,21,41,61,81,101] # e.g. 1 for 0, 11 for 0.1, 51 for 0.5, 101 for 1.

# Flag to specify if strategy outputs should be computed and plotted
compute_optim_strategy_flag = false
plot_optim_strategy_flag = false

# Flag to specify if cost outputs should be produced
plot_costs_flag = false

# Specify if infection & vaccination outputs should be plotted and what
# percentile values + exceedence of threshold infection should be computed
plot_infections_vaccs_flag = false
compute_infections_vacc_flag = false # If false, values are loaded from an existing MAT file, selected based on config_ID
prctiles_infections_vaccs = [2.5,50,97.5]
threshold_infections_vals = [10,25,50,100,250,500,1000]

# Call the function to produce the ternary plot set
percentage_threshold_inf_exceeded_array_animallevel_popnpersp,
  percentage_threshold_inf_exceeded_array_animallevel_indivpersp,
  percentage_infected_array_animallevel_popnpersp,
  percentage_infected_array_animallevel_indivpersp,
  percentage_vacc_array_animallevel_popnpersp,
  percentage_vacc_array_animallevel_indivpersp = OptimThresholdVsCostRatio_SellkeRuns_GenerateTernaryPlot(save_fig_flag,
                                                            config_ID,
                                                            scen_IDs,
                                                            batch_ID_offset,
                                                            rel_cost_of_vacc_idxs,
                                                            compute_optim_strategy_flag,
                                                            plot_optim_strategy_flag,
                                                            plot_costs_flag,
                                                            plot_infections_vaccs_flag,
                                                            compute_infections_vacc_flag,
                                                            prctiles_infections_vaccs,
                                                            threshold_infections_vals)
# Check if should save outputs to MAT file
if compute_infections_vacc_flag == true

     # Set save filename based on configuration
    if config_ID == "cumbria"
        save_file_MAT = "OptimBehaviour_InfVacc_PlotData/cumbria_inf_vacc_ternary_plot_data_julia.mat"
    elseif config_ID == "cumbria_alt"
        save_file_MAT = "OptimBehaviour_InfVacc_PlotData/cumbria_alt_inf_vacc_ternary_plot_data_julia.mat"
    elseif config_ID == "devon"
        save_file_MAT = "OptimBehaviour_InfVacc_PlotData/devon_inf_vacc_ternary_plot_data_julia.mat"
    elseif config_ID == "devon_alt"
        save_file_MAT = "OptimBehaviour_InfVacc_PlotData/devon_alt_inf_vacc_ternary_plot_data_julia.mat"
    else
        error("Invalid config_ID provided.")
    end

    matwrite(save_file_MAT, Dict(
            "percentage_threshold_inf_exceeded_array_animallevel_popnpersp" => percentage_threshold_inf_exceeded_array_animallevel_popnpersp,
            "percentage_threshold_inf_exceeded_array_animallevel_indivpersp" => percentage_threshold_inf_exceeded_array_animallevel_indivpersp,
            "percentage_infected_array_animallevel_popnpersp" => percentage_infected_array_animallevel_popnpersp,
            "percentage_infected_array_animallevel_indivpersp" => percentage_infected_array_animallevel_indivpersp,
            "percentage_vacc_array_animallevel_popnpersp" => percentage_vacc_array_animallevel_popnpersp,
            "percentage_vacc_array_animallevel_indivpersp" => percentage_vacc_array_animallevel_indivpersp
            ); compress = true)
end
