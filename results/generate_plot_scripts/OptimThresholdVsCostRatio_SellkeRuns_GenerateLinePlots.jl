#=
Purpose:
Produce line plots of optimal decision point against cost ratio. Also generates
infection and vaccination summary statistic plots.
=#

#-------------------------------------------------------------------------------
# SUPPORTING FUNCTIONS
#-------------------------------------------------------------------------------
function modify_value!(data_array,original_val,plot_val)

    # For values of data_array matching original_val, replace with plot_val
    data_array[data_array.==original_val] .= plot_val

    return nothing
end


#-------------------------------------------------------------------------------
# MAIN FUNCTION
#-------------------------------------------------------------------------------
"""
    OptimThresholdVsCostRatio_SellkeRuns_GenerateLinePlots(save_fig_flag::Bool,
                                                            config_ID::String)

Produce line plots of optimal decision point against cost ratio. Also generates
infection and vaccination summary statistic plots.

Inputs:
- `save_fig_flag::Bool`: Flag determining whether generated figures are saved to file or not.
- `config_ID::String`: The demography, location and epidemiological configuration in use.

Outputs: None \n
Location: OptimThresholdVsCostRatio\\_SellkeRuns\\_GenerateLinePlots.jl
"""
function OptimThresholdVsCostRatio_SellkeRuns_GenerateLinePlots(save_fig_flag::Bool,
                                                                config_ID::String)

      #-------------------------------------------------------------------------
      # Specify variables dependent on configuration in use
      # Input filename specifier, save filename ID, plot titles & legend inclusion
      #-------------------------------------------------------------------------
      if config_ID == "cumbria"
         input_filename_suffix = "Cumbria_vacc_distance_risk_measure_scenID921"
         save_filename_suffix = "OptimBehaviour_SellkeRuns_plots_julia/cumbria_scenID921_"
         county_specific_plot_title = "Cumbria"
         include_legend_flag = false
      elseif config_ID == "cumbria_alt"
         input_filename_suffix = "Cumbria_alt_pathogen_vacc_distance_risk_measure_scenID1221"
         save_filename_suffix = "OptimBehaviour_SellkeRuns_plots_julia/cumbria_alt_scenID1221_"
         county_specific_plot_title = "Cumbria"
         include_legend_flag = false
      elseif config_ID == "devon"
         input_filename_suffix = "Devon_vacc_distance_risk_measure_scenID921"
         save_filename_suffix = "OptimBehaviour_SellkeRuns_plots_julia/devon_scenID921_"
         county_specific_plot_title = "Devon"
         include_legend_flag = true
      elseif config_ID == "devon_alt"
         input_filename_suffix = "Devon_alt_pathogen_vacc_distance_risk_measure_scenID1221"
         save_filename_suffix = "OptimBehaviour_SellkeRuns_plots_julia/devon_alt_scenID1221_"
         county_specific_plot_title = "Devon"
         include_legend_flag = true
      else
         error("Invalid config_ID provided.")
      end

      OptimThresholdVals_PremLevel_PopnPersp_FileName = string("OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_PremLevel_PopnPersp_", input_filename_suffix, ".txt")
      OptimThresholdVals_AnimalLevel_PopnPersp_FileName = string("OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_PopnPersp_", input_filename_suffix,".txt")

      OptimThresholdVals_PremLevel_IndivPersp_FileName = string("OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_PremLevel_IndivPersp_", input_filename_suffix, ".txt")
      OptimThresholdVals_AnimalLevel_IndivPersp_FileName = string("OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_IndivPersp_", input_filename_suffix,".txt")

      #--------------------------------------------------------------------------
      # Load summary statistic data
      #--------------------------------------------------------------------------
      OptimThresholdVals_PremLevel_PopnPersp = readdlm(OptimThresholdVals_PremLevel_PopnPersp_FileName,'\t')
      OptimThresholdVals_AnimalLevel_PopnPersp = readdlm(OptimThresholdVals_AnimalLevel_PopnPersp_FileName,'\t')

      OptimThresholdVals_PremLevel_IndivPersp = readdlm(OptimThresholdVals_PremLevel_IndivPersp_FileName,'\t')
      OptimThresholdVals_AnimalLevel_IndivPersp = readdlm(OptimThresholdVals_AnimalLevel_IndivPersp_FileName,'\t')

      #----------------------------------------------------------------------
      # Load infection and vaccination (under optimal strategy) summary statistic data
      #----------------------------------------------------------------------
      if config_ID == "cumbria"
         infection_vacc_data = matread("OptimBehaviour_InfVacc_PlotData/cumbria_inf_vacc_ternary_plot_data.mat")
      elseif config_ID == "cumbria_alt"
         infection_vacc_data = matread("OptimBehaviour_InfVacc_PlotData/cumbria_alt_inf_vacc_ternary_plot_data.mat")
      elseif config_ID == "devon"
         infection_vacc_data = matread("OptimBehaviour_InfVacc_PlotData/devon_inf_vacc_ternary_plot_data.mat")
      elseif config_ID == "devon_alt"
         infection_vacc_data = matread("OptimBehaviour_InfVacc_PlotData/devon_alt_inf_vacc_ternary_plot_data.mat")
      else
          error("Invalid config_ID provided.")
      end

      # Assign relevant variables from data structure to plotting variables
      # Column 21 for scenario 21 (all reactionary scenario)
          # Infection, slice to pick out the threshold count to match or exceed: threshold_infections_vals = [10 25 50 100 250 500 1000];
      infection_exceed_thereshold_popnpersp = infection_vacc_data["percentage_threshold_inf_exceeded_array_animallevel_popnpersp"][:,21,2]
      infection_exceed_thereshold_indivpersp = infection_vacc_data["percentage_threshold_inf_exceeded_array_animallevel_indivpersp"][:,21,2]

          # Percentage of premises infected, slice to pick out the percentile: [2.5 50 97.5] e.g. Index 2 for median.
      infection_data_popnpersp_median = infection_vacc_data["percentage_infected_array_animallevel_popnpersp"][:,21,2]
      infection_data_indivpersp_median = infection_vacc_data["percentage_infected_array_animallevel_indivpersp"][:,21,2]

      infection_data_popnpersp_LB = infection_vacc_data["percentage_infected_array_animallevel_popnpersp"][:,21,1]
      infection_data_indivpersp_LB = infection_vacc_data["percentage_infected_array_animallevel_indivpersp"][:,21,1]

      infection_data_popnpersp_UB = infection_vacc_data["percentage_infected_array_animallevel_popnpersp"][:,21,3]
      infection_data_indivpersp_UB = infection_vacc_data["percentage_infected_array_animallevel_indivpersp"][:,21,3]

          # Percentage of premises vaccinated, slice to pick out the percentile: [2.5 50 97.5] e.g. Index 2 for median.
      vacc_data_popnpersp_median = infection_vacc_data["percentage_vacc_array_animallevel_popnpersp"][:,21,2]
      vacc_data_indivpersp_median = infection_vacc_data["percentage_vacc_array_animallevel_indivpersp"][:,21,2]

      vacc_data_popnpersp_LB = infection_vacc_data["percentage_vacc_array_animallevel_popnpersp"][:,21,1]
      vacc_data_indivpersp_LB = infection_vacc_data["percentage_vacc_array_animallevel_indivpersp"][:,21,1]

      vacc_data_popnpersp_UB = infection_vacc_data["percentage_vacc_array_animallevel_popnpersp"][:,21,3]
      vacc_data_indivpersp_UB = infection_vacc_data["percentage_vacc_array_animallevel_indivpersp"][:,21,3]

      #--------------------------------------------------------------------------
      # Declare vaccine to infection cost ratio values that were tested
      #--------------------------------------------------------------------------
      VaccToInfCostRatio = 0:0.01:1  # We set C_I = 1

      #--------------------------------------------------------------------------
      # Set up plotting variables
      #--------------------------------------------------------------------------
      # Axes labels
      yaxis_name = L"\textrm{Notified \ premises \ distance \ threshold \ (km)}"
      yaxis_name_infection = L"\textrm{Replicates \ with \ 25+ \ premises \ infected \ (\%)}"
      yaxis_name_infection_with_PI = L"\textrm{Premises \ infected \ (\%)}"
      yaxis_name_vacc = L"\textrm{Premises \ vaccinated \ (\%)}"
      # yaxis_name = "Risk threshold: Notified premises distance (km)"
      # yaxis_name_infection = "Replicates with 25+ premises infected (%)"
      # yaxis_name_infection_with_PI = "Premises infected (%)"
      # yaxis_name_vacc = "Premises vaccinated (%)"

      xaxis_name = L"\textrm{Relative \ cost \ of \ vaccination,} \mathrm{C_V}"

      # y-axis limits & tick labels - optimal strategy plots
      yaxis_lims = (0,10)
      yaxis_tick_pos = 0:1:10
      yaxis_tick_labels = ["0","1","2","3","4","5","6","7","8","9","10"]

      # y-axis limits & tick labels - infection plots (median profile only)
      yaxis_lims_infection = (0,100)
      yaxis_tick_pos_infection = 0:10:100
      yaxis_tick_labels_infection = 0:10:100

      # y-axis limits & tick labels - infection plots (with prediction intervals)
      yaxis_lims_infection_with_PI = (0,42)
      yaxis_tick_pos_infection_with_PI = 0:4:40
      yaxis_tick_labels_infection_with_PI = 0:4:40

      # y-axis limits & tick labels - vaccination plots (median profile only)
      yaxis_lims_vacc = (0,10)
      yaxis_tick_pos_vacc = 0:1:10
      yaxis_tick_labels_vacc = 0:1:10

      # y-axis limits & tick labels - vaccination plots (with prediction intervals)
      yaxis_lims_vacc_with_PI = (0,20)
      yaxis_tick_pos_vacc_with_PI = 0:2:20
      yaxis_tick_labels_vacc_with_PI = 0:2:20

      # x-axis tick labels
      xaxis_lims = (0,1)
      xaxis_tick_pos = 0:0.1:1
      xaxis_tick_labels = ["0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"]

      # Set plot face transparency
      facealpha_val = 0.2

      # Fontsize
      label_fontsize = 15
      tick_fontsize = 13

      # Set figure width and height
      fig_width = 550
      fig_height = 450

      #--------------------------------------------------------------------------
      # Assign select data to variables
      #--------------------------------------------------------------------------

      # Note, summary statistics output are the mean and the following quantiles:
      # [mean,0.025,0.5,0.975,0,1,0.25,0.75,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9]

      # Pick out median values
      OptimThresholdVals_PremLevel_PopnPersp_Median = OptimThresholdVals_PremLevel_PopnPersp[:,3]
      OptimThresholdVals_AnimalLevel_PopnPersp_Median = OptimThresholdVals_AnimalLevel_PopnPersp[:,3]

      OptimThresholdVals_PremLevel_IndivPersp_Median = OptimThresholdVals_PremLevel_IndivPersp[:,3]
      OptimThresholdVals_AnimalLevel_IndivPersp_Median = OptimThresholdVals_AnimalLevel_IndivPersp[:,3]

      # #Get lower bound and upper bounds
      # OptimThresholdVals_PremLevel_PopnPersp_95PI_LB = OptimThresholdVals_PremLevel_PopnPersp[:,2]
      # OptimThresholdVals_PremLevel_PopnPersp_95PI_UB = OptimThresholdVals_PremLevel_PopnPersp[:,4]
      #
      # OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB = OptimThresholdVals_AnimalLevel_PopnPersp[:,2]
      # OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB = OptimThresholdVals_AnimalLevel_PopnPersp[:,4]
      #
      # OptimThresholdVals_PremLevel_IndivPersp_95PI_LB = OptimThresholdVals_PremLevel_IndivPersp[:,2]
      # OptimThresholdVals_PremLevel_IndivPersp_95PI_UB = OptimThresholdVals_PremLevel_IndivPersp[:,4]
      #
      # OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB = OptimThresholdVals_AnimalLevel_IndivPersp[:,2]
      # OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB = OptimThresholdVals_AnimalLevel_IndivPersp[:,4]


      #--------------------------------------------------------------------------
      # Construct optimal control line profiles (displays costs applied at animal level only)
      #--------------------------------------------------------------------------
      # Construct line plots
      p1 = plot(VaccToInfCostRatio,
                  OptimThresholdVals_AnimalLevel_PopnPersp_Median,
                  linestyle = :solid,
                  color = RGB(0,0,0.8),
                  linewidth=1.5,
                  framestyle = :box,
                  fontfamily = "Computer Modern",
                  label = "Population perspective"
                  )

      plot!(VaccToInfCostRatio,
            OptimThresholdVals_AnimalLevel_IndivPersp_Median,
            linestyle = :dash,
            color = RGB(0.8,0,0),
            linewidth=1.5,
            label = "Individual perspective")

     # Add title
    plot!(title = county_specific_plot_title,
            titlefont = "Helvetica Bold")

     # Alter x-axis properties
     plot!(xlabel = xaxis_name,
            xticks = (xaxis_tick_pos, xaxis_tick_labels),
            xlims = xaxis_lims,
            xguidefontsize = label_fontsize,
            xtickfontsize = tick_fontsize,
            xformatter = :plain)

     # Alter y-axis properties
     plot!(ylabel = yaxis_name,
             yticks = (yaxis_tick_pos, yaxis_tick_labels),
             ylims = yaxis_lims,
             yguidefontsize = label_fontsize,
             ytickfontsize = tick_fontsize,
             yformatter = :plain)

      # Add legend (if applicable)
      plot!(legend = include_legend_flag,
            legendfontsize = tick_fontsize)

      # If appropriate, specify the overall plot size
      plot!(size=(fig_width, fig_height))

      # Save figure
      if save_fig_flag == true
          save_filename = string(save_filename_suffix,"optimise_median_cost_animallevel_sellke_julia.pdf")
          savefig(p1, save_filename)
      end


      #----------------------------------------------------------------------
      # Construct infection line profiles under optimal control strategy with
      # prediction intervals
      # (displays costs applied at animal level only)
      #----------------------------------------------------------------------

      # To shade region between two curves, use fillrange
      # fillrange: plot(x, l, fillrange = u, fillalpha = 0.35, c = 1, label = "Confidence band")

      # Construct line for median population perspective
      p2 = plot(VaccToInfCostRatio,
                  infection_data_popnpersp_median,
                  linestyle = :solid,
                  color = RGB(0,0,0.8),
                  linewidth=1.5,
                  framestyle = :box,
                  fontfamily = "Computer Modern",
                  label = ""
                  )

      # Construct line for median individual perspective
      plot!(VaccToInfCostRatio,
            infection_data_indivpersp_median,
            linestyle = :dash,
            color = RGB(0.8,0,0),
            linewidth=1.5,
            label = "")

      # Construct filled region - population perspective
      # syntax: plot(x, l, fillrange = u, fillalpha = 0.35, fillcolor = 1, label = "Confidence band")
      plot!(VaccToInfCostRatio, infection_data_popnpersp_LB,
                fillrange = infection_data_popnpersp_UB,
                fillalpha = facealpha_val,
                linealpha = facealpha_val,
                linewidth = 0,
                c = RGB(0,0,0.8),
                label = "",
                )

      # Construct filled region - individual perspective
      # syntax: plot(x, l, fillrange = u, fillalpha = 0.35, fillcolor = 1, label = "Confidence band")
      plot!(VaccToInfCostRatio, infection_data_indivpersp_LB,
              fillrange = infection_data_indivpersp_UB,
              fillalpha = facealpha_val,
              linealpha = facealpha_val,
              linewidth = 0,
              c = RGB(0.8,0,0),
              label = "",
              )

     # Add title
     plot!(title = county_specific_plot_title,
            titlefont = "Helvetica Bold")

     # Alter x-axis properties
     plot!(xlabel = xaxis_name,
            xticks = (xaxis_tick_pos, xaxis_tick_labels),
            xlims = xaxis_lims,
            xguidefontsize = label_fontsize,
            xtickfontsize = tick_fontsize,
            xformatter = :plain)

     # Alter y-axis properties
     plot!(ylabel = yaxis_name_infection_with_PI,
             yticks = (yaxis_tick_pos_infection_with_PI, yaxis_tick_labels_infection_with_PI),
             ylims = yaxis_lims_infection_with_PI,
             yguidefontsize = label_fontsize,
             ytickfontsize = tick_fontsize,
             yformatter = :plain)

      # If appropriate, specify the overall plot size
      plot!(size=(fig_width, fig_height))

      # Save figure
      if save_fig_flag == true
          save_filename = string(save_filename_suffix,"all_reactionary_infection_profile_with_PI_sellke_julia.pdf")
          savefig(p2, save_filename)
      end

      #----------------------------------------------------------------------
      # Construct vaccination line profiles under optimal control strategy with
      # prediction intervals
      # (displays costs applied at animal level only)
      #----------------------------------------------------------------------

      # To shade region between two curves, use fillrange
      # syntax: plot([x x], fillrange = [l u], fillalpha = 0.35, fillcolor = 1, label = "Confidence band")

      # Construct line for median population perspective
      p3 = plot(VaccToInfCostRatio,
                  vacc_data_popnpersp_median,
                  linestyle = :solid,
                  color = RGB(0,0,0.8),
                  linewidth = 1.5,
                  framestyle = :box,
                  fontfamily = "Computer Modern",
                  label = ""
                  )

      # Construct line for median individual perspective
      plot!(VaccToInfCostRatio,
            vacc_data_indivpersp_median,
            linestyle = :dash,
            color = RGB(0.8,0,0),
            linewidth = 1.5,
            label = "")

      # Construct filled region - population perspective
      # syntax: plot(x, l, fillrange = u, fillalpha = 0.35, fillcolor = 1, label = "Confidence band")
      plot!(VaccToInfCostRatio, vacc_data_popnpersp_LB,
                fillrange = vacc_data_popnpersp_UB,
                fillalpha = facealpha_val,
                linealpha = facealpha_val,
                c = RGB(0,0,0.8),
                label = "",
                )

      # Construct filled region - individual perspective
      # syntax: plot(x, l, fillrange = u, fillalpha = 0.35, fillcolor = 1, label = "Confidence band")
      plot!(VaccToInfCostRatio, vacc_data_indivpersp_LB,
              fillrange = vacc_data_indivpersp_UB,
              fillalpha = facealpha_val,
              linealpha = facealpha_val,
              c = RGB(0.8,0,0),
              label = "",
              )

      # Add title
      plot!(title = county_specific_plot_title,
            titlefont = "Helvetica Bold")

      # Alter x-axis properties
      plot!(xlabel = xaxis_name,
            xticks = (xaxis_tick_pos, xaxis_tick_labels),
            xlims = xaxis_lims,
            xguidefontsize = label_fontsize,
            xtickfontsize = tick_fontsize,
            xformatter = :plain)

      # Alter y-axis properties
      plot!(ylabel = yaxis_name_vacc,
             yticks = (yaxis_tick_pos_vacc_with_PI, yaxis_tick_labels_vacc_with_PI),
             ylims = yaxis_lims_vacc_with_PI,
             yguidefontsize = label_fontsize,
             ytickfontsize = tick_fontsize,
             yformatter = :plain)

      # If appropriate, specify the overall plot size
      plot!(size=(fig_width, fig_height))

      # Save figure
      if save_fig_flag == true
          save_filename = string(save_filename_suffix,"all_reactionary_vacc_profile_with_PI_sellke_julia.pdf")
          savefig(p3, save_filename)
      end

      return nothing
end

# #--------------------------------------------------------------------------
# # Construct plot, optimal strategy based on 2.5th percentile
# #--------------------------------------------------------------------------
#
# # Construct line plots
# p3 = plot(VaccToInfCostRatio,
#             OptimThresholdVals_PremLevel_PopnPersp_95PI_LB,
#             color = RGB(0,0,0.8),
#             linewidth=1.5,
#             framestyle = :box,
#             title = "Premises-level costs",
#             legend = false)
#
# plot!(VaccToInfCostRatio,
#       OptimThresholdVals_PremLevel_IndivPersp_95PI_LB,
#       linestyle = :dash,
#       color = RGB(0.8,0,0),
#       linewidth=1.5)
#
# p4 = plot(VaccToInfCostRatio,
#             OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB,
#             linestyle = :solid,
#             color = RGB(0,0,0.8),
#             linewidth=1.5,
#             framestyle = :box,
#             title = "Animal-level costs",
#             legend = false,
#             label = "Population perspective")
#
# plot!(VaccToInfCostRatio,
#       OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB,
#       linestyle = :dash,
#       color = RGB(0.8,0,0),
#       linewidth=1.5,
#       label = "Individual perspective")
#
# # Construct plot layout
# p_LB = plot(p3, p4,
#              layout = (1, 2),
#              xformatter = :plain,
#              yformatter = :plain,
#              xlabel = xaxis_name,
#              ylabel = yaxis_name,
#              yticks = (yaxis_tick_pos, yaxis_tick_labels),
#              ylims = yaxis_lims)
#
# # If appropriate, specify the overall plot size
# # plot!(size=(750, 450))
#
# # Save figure
# savefig(p_LB, "OptimBehaviour_Plots/CumbriaExamplePlot_optimise_LB_percentile_julia.pdf")
#
# #--------------------------------------------------------------------------
# # Construct plot, optimal strategy based on 97.5th percentile
# #--------------------------------------------------------------------------
#
# # Construct line plots
# p5 = plot(VaccToInfCostRatio,
#             OptimThresholdVals_PremLevel_PopnPersp_95PI_UB,
#             color = RGB(0,0,0.8),
#             linewidth=1.5,
#             framestyle = :box,
#             title = "Premises-level costs",
#             legend = false)
#
# plot!(VaccToInfCostRatio,
#       OptimThresholdVals_PremLevel_IndivPersp_95PI_UB,
#       linestyle = :dash,
#       color = RGB(0.8,0,0),
#       linewidth=1.5)
#
# p6 = plot(VaccToInfCostRatio,
#             OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB,
#             linestyle = :solid,
#             color = RGB(0,0,0.8),
#             linewidth=1.5,
#             framestyle = :box,
#             title = "Animal-level costs",
#             legend = false,
#             label = "Population perspective")
#
# plot!(VaccToInfCostRatio,
#       OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB,
#       linestyle = :dash,
#       color = RGB(0.8,0,0),
#       linewidth=1.5,
#       label = "Individual perspective")
#
# # Construct plot layout
# p_UB = plot(p5, p6,
#              layout = (1, 2),
#              xformatter = :plain,
#              yformatter = :plain,
#              xlabel = xaxis_name,
#              ylabel = yaxis_name,
#              yticks = (yaxis_tick_pos, yaxis_tick_labels),
#              ylims = yaxis_lims)
#
# # If appropriate, specify the overall plot size
# # plot!(size=(1100, 450))
#
# # Save figure
# savefig(p_UB, "OptimBehaviour_Plots/CumbriaExamplePlot_optimise_UB_percentile_julia.pdf")
