#=
Purpose:
Produce outputs of optimal strategy (out of a collection of strategies)
based on three different groups
 (i) Initial control adopters (enact control from outset of simualation)
 (ii) Reactive responders (implement control based on epidemic situation)
 (iii) Never implement control

Plot produced for a single relative cost of vaccination
Outbreak simulations used the Sellke construction.

Also can produce ternary plots for:
 - cost under optimal strategy
 - estimated infections (for specified percentiles) under optimal strategy
 - estimated vaccinations (for specified percentiles) under optimal strategy
=#

"""
    OptimThresholdVsCostRatio_SellkeRuns_GenerateTernaryPlot(save_figure_flag::Bool,
                                                                    config_ID::String,
                                                                    scen_IDs::Array{Int64,1},
                                                                    batch_ID_offset::Int64,
                                                                    rel_cost_of_vacc_idxs::Array{Int64,1},
                                                                    compute_optim_strategy_flag::Bool,
                                                                    plot_optim_strategy_flag::Bool,
                                                                    plot_costs_flag::Bool,
                                                                    plot_infections_vaccs_flag::Bool,
                                                                    compute_infections_vacc_flag::Bool,
                                                                    prctiles_infections_vaccs::Array{Float64,1},
                                                                    threshold_infections_vals::Array{Int64,1})

Produce ternary plots of optimal decision point against cost ratio.
Also generates infection and vaccination summary ternary plots.

Inputs:
- `save_fig_flag::Bool`: Flag determining whether generated figures are saved to file or not.
- `config_ID::String`: The demography, location and epidemiological configuration in use.
- `scen_IDs::Array{Int64,1}`: List of scenarios to be plotted
- `batch_ID_offset::Int64`: Added to the batch_ID (for accessing model output files)
- `rel_cost_of_vacc_idxs::Array{Int64,1}`: Relative cost of vaccination to generate plots for
- `compute_optim_strategy_flag::Bool`: Flag to specify if strategy outputs should be computed
- `plot_optim_strategy_flag::Bool`: Flag to specify if strategy outputs should be plotted
- `plot_costs_flag::Bool`: Flag to specify if cost outputs should be produced
- `plot_infections_vaccs_flag::Bool`: Flag to specify if infection & vaccination outputs should be plotted
- `compute_infections_vacc_flag::Bool`: Flag to specify if infection & vaccination outputs should be computed
- `prctiles_infections_vaccs::Array{Float64,1}`: percentile values to be computed
- `threshold_infections_vals::Array{Int64,1}`: pexceedence of threshold infection to be computed

Outputs: For different perspectives, percentages for replicates exceeding infection threshold, premises infected & vaccinated under the optimal strategy
- `percentage_threshold_inf_exceeded_array_animallevel_popnpersp::Array{Float64,3}`
- `percentage_threshold_inf_exceeded_array_animallevel_indivpersp::Array{Float64,3}`
- `percentage_infected_array_animallevel_popnpersp::Array{Float64,3}`
- `percentage_infected_array_animallevel_indivpersp::Array{Float64,3}`
- `percentage_vacc_array_animallevel_popnpersp::Array{Float64,3}`
- `percentage_vacc_array_animallevel_indivpersp::Array{Float64,3}`

Location: OptimThresholdVsCostRatio\\_SellkeRuns\\_GenerateTernaryPlots.jl
"""
function OptimThresholdVsCostRatio_SellkeRuns_GenerateTernaryPlot(save_figure_flag::Bool,
                                                                    config_ID::String,
                                                                    scen_IDs::Array{Int64,1},
                                                                    batch_ID_offset::Int64,
                                                                    rel_cost_of_vacc_idxs::Array{Int64,1},
                                                                    compute_optim_strategy_flag::Bool,
                                                                    plot_optim_strategy_flag::Bool,
                                                                    plot_costs_flag::Bool,
                                                                    plot_infections_vaccs_flag::Bool,
                                                                    compute_infections_vacc_flag::Bool,
                                                                    prctiles_infections_vaccs::Array{Float64,1},
                                                                    threshold_infections_vals::Array{Int64,1})


    #--------------------------------------------------------------------------
    # Set ternary plot variables
    #--------------------------------------------------------------------------

    # Specify number of control strategies considered
    n_risk_measures = 11 # 0km, 1km, . , 9km, 10km distance based measures

    # Axes labels
    ternary_axes_labels = ["Precautionary (%)","Reactionary (%)","Non-vaccinators (%)"]

    # Font sizes
    label_fontsize = 12
    tick_label_fontsize = 10

    # Colourbar
    optim_strat_cbar_range = (-0.5,n_risk_measures-0.5)
    optim_strat_colourbar_tick_pos = collect(0:1.:(n_risk_measures-1))
    optim_strat_colourbar_tick_labels = string.(collect(0:1:10))
    optim_strat_colourbar_label = "Risk threshold: Notified premises distance (km)"

    # Setup colour map options for optimal strategy plots
    CT_red = cgrad(:Reds_9, n_risk_measures, categorical = true)

    my_cmap = [RGB(0.5,0.5,0.5);
                         CT_red[end:-1:2]] # Have first entry as grey
    optim_strat_chosen_cmap = palette(my_cmap, n_risk_measures) # Set uo colourmap as a palette to give discrete colourbar (rather than a gradient)

    # # Setup colourbar options for cost under optimal strategy plots
    # CT_blue_cost_plots = cgrad(:Blues_9)
    # strat_cost_chosen_cmap = CT_blue_cost_plots
    # strat_cost_colourbar_label_popn_persp = "Population perspective cost (using optimal strategy)"
    # strat_cost_colourbar_label_indiv_persp = "Individual perspective cost (using optimal strategy)"
    # strat_cost_cbar_range = [0 1500000]
    # strat_cost_colourbar_tick_pos = [0 250000 500000 750000 1000000 1250000 1500000]
    # strat_cost_colourbar_tick_labels = 0:250000:1500000

    # Setup colourbar options for exceeding threshold infection value under optimal strategy plots
    strat_exceed_threshold_infections_chosen_cmap = cgrad(:Oranges_9)
    strat_exceed_threshold_infections_colourbar_label = "Replicates with 25+ premises infected (%)"
    strat_exceed_threshold_infections_cbar_range = (0.,100.)
    strat_exceed_threshold_infections_colourbar_tick_pos = collect(0:10.:100)
    strat_exceed_threshold_infections_colourbar_tick_labels = string.(collect(0:10:100))

     # # Setup colourbar options for percentage infections under optimal strategy plots
     # strat_infections_chosen_cmap = cgrad(:Oranges_9)
     # strat_infections_colourbar_label = "Premises infected (%)"
     # strat_infections_cbar_range = [0 10]
     # strat_infections_colourbar_tick_pos = 0:1:10
     # strat_infections_colourbar_tick_labels = ["0","1","2","3","4","5","6","7","8","9","10+"]

    # Setup colourbar options for vaccinations under optimal strategy plots
    strat_vaccs_chosen_cmap = cgrad(:Blues_9)
    strat_vaccs_colourbar_label = "Premises vaccinated (%)"
    strat_vaccs_cbar_range = (0.,100.)
    strat_vaccs_colourbar_tick_pos = [0.,10,20,30,40,50,60,70,80,90,100]
    strat_vaccs_colourbar_tick_labels = string.([0,10,20,30,40,50,60,70,80,90,100])

    # Specify amount of infection related thresholds in use
    n_prctiles_infections_vaccs = length(prctiles_infections_vaccs)
    n_threshold_infections_vals = length(threshold_infections_vals)

    #--------------------------------------------------------------------------
    # Load optimum strategy data
    #--------------------------------------------------------------------------

    # Get number of scenarios
    n_scens = length(scen_IDs)

    # Set aggregated epidemiological output file directory to be accessed
    # and number of premises in use
    if (config_ID == "cumbria") || (config_ID == "cumbria_alt")
       epi_outputs_agg_directory = "../GBCountyModelSimnOutputs/Cumbria_EpiOutputs_Aggregated/"
       prem_num = 3784
    elseif (config_ID == "devon") || (config_ID == "devon_alt")
       epi_outputs_agg_directory = "../GBCountyModelSimnOutputs/Devon_EpiOutputs_Aggregated/"
       prem_num = 5343
    else
       error("Invalid config_ID provided.")
    end

    # Initialise filename vector
    filename_suffix_vec = Array{String,1}(undef,n_scens)

    # Populate filename vector
    for scen_itr = 1:n_scens

      # Get scenario ID for current iteration
      current_itr_scen_ID = scen_IDs[scen_itr]

      # Assign filename to cel
      if config_ID == "cumbria"
          filename_suffix_vec[scen_itr] = "Cumbria_vacc_distance_risk_measure_scenID$current_itr_scen_ID"
       elseif config_ID == "cumbria_alt"
          filename_suffix_vec[scen_itr] = "Cumbria_alt_pathogen_vacc_distance_risk_measure_scenID$current_itr_scen_ID"
       elseif config_ID == "devon"
          filename_suffix_vec[scen_itr] = "Devon_vacc_distance_risk_measure_scenID$current_itr_scen_ID"
       elseif config_ID == "devon_alt"
          filename_suffix_vec[scen_itr] = "Devon_alt_pathogen_vacc_distance_risk_measure_scenID$current_itr_scen_ID"
      else
          error("Invalid config_ID provided.")
      end
    end

    # Assign number of vaccine group combinations tested to variable
    n_vacc_uptake_strats = length(filename_suffix_vec)

    #--------------------------------------------------------------------------
    # Construct ternary plot position data
    #--------------------------------------------------------------------------

    # Position array format: Column 1: Initial adopters fraction
    #                        Column 2: Reactioners fraction
    #                        Column 3: Non-vaccinators fraction
    ternay_plot_pos_array = zeros(Float64,n_scens,3)

    # Load group composition per scenario from data file
    control_param_input_args_array = readdlm("../../../src/SpatialSimns/ToyModel/parameter_combination_files/generic_model_control_param_val_array.txt")
        # Three column array.
        # Column 1: Proportion in "Initial adopters" group
        # Column 2: Proportion inf "Reactioners" group
        # Column 3: Risk threshold measure (NOT USED HERE)

    # Assign "Initial adopters" & "Reactioners" group columns to array
    # Take every "n_risk_measures"th row from control_param_input_args_array,
    # as that corresponds to moving to the next scenario chunk that used a
    # different group composition
    ternay_plot_pos_array[:,1:2] = control_param_input_args_array[1:n_risk_measures:end,1:2]

    # Populate "No control" group using data on other groups
    ternay_plot_pos_array[:,3] = ones(Float64,n_scens) - sum(ternay_plot_pos_array[:,1:2],dims=2)

    # Error check
    if sum((sum(ternay_plot_pos_array,dims=2) .!= 1)) > 0
        error("Rows in ternay_plot_pos_array do not sum to 1 as expected.")
    end

    #--------------------------------------------------------------------------
    # Construct optimal strategy arrays for all relative vaccine costs
    # Assign the cost under the optimal strategy to array
    #--------------------------------------------------------------------------

    # Specify the number of relative vaccination costs tested
    n_rel_vacc_costs = 101

    # Specify the percentile to be used in plots
    # Slice idx: 1 - mean; 2 - 2.5th percentile; 3 - median; 4 - 97.5th
    # percentile; 5 - minimum; 6 - maximum.
    percentile_slice_idx = 3

    # Initialise storage arrays for the optimal strategy
    optim_control_array_animallevel_popnpersp = zeros(Int64,n_rel_vacc_costs,n_vacc_uptake_strats)
    optim_control_array_animallevel_indivpersp = zeros(Int64,n_rel_vacc_costs,n_vacc_uptake_strats)

    # Initialise storage arrays to record the percentage of replicates that
    # exceeded a specified number of premises infected under the strategy identified as optimal
    percentage_threshold_inf_exceeded_array_animallevel_popnpersp = zeros(Float64,n_rel_vacc_costs,n_vacc_uptake_strats,n_threshold_infections_vals)
    percentage_threshold_inf_exceeded_array_animallevel_indivpersp = zeros(Float64,n_rel_vacc_costs,n_vacc_uptake_strats,n_threshold_infections_vals)

    # Initialise storage arrays to record the percentage of premises infected
    # under the strategy identified as optimal
    percentage_infected_array_animallevel_popnpersp = zeros(Float64,n_rel_vacc_costs,n_vacc_uptake_strats,n_prctiles_infections_vaccs)
    percentage_infected_array_animallevel_indivpersp = zeros(Float64,n_rel_vacc_costs,n_vacc_uptake_strats,n_prctiles_infections_vaccs)

    # Initialise storage arrays to record the percentage of premises vaccinated
    # under the strategy identified as optimal
    percentage_vacc_array_animallevel_popnpersp = zeros(Float64,n_rel_vacc_costs,n_vacc_uptake_strats,n_prctiles_infections_vaccs)
    percentage_vacc_array_animallevel_indivpersp = zeros(Float64,n_rel_vacc_costs,n_vacc_uptake_strats,n_prctiles_infections_vaccs)

    # Populate arrays, iterate over each scenario
    for vacc_strat_itr = 1:n_vacc_uptake_strats

        #=
        Get the optimal strategy data
        =#
        if compute_optim_strategy_flag == true
            # Set relevant data file names
            OptimThresholdVals_AnimalLevel_PopnPersp_FileName = "OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_PopnPersp_$(filename_suffix_vec[vacc_strat_itr]).txt"
            OptimThresholdVals_AnimalLevel_IndivPersp_FileName = "OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_IndivPersp_$(filename_suffix_vec[vacc_strat_itr]).txt"

            # Load optimal strategy data
            OptimThresholdVals_AnimalLevel_PopnPersp = readdlm(OptimThresholdVals_AnimalLevel_PopnPersp_FileName)
            OptimThresholdVals_AnimalLevel_IndivPersp = readdlm(OptimThresholdVals_AnimalLevel_IndivPersp_FileName)

            # Pick out relevant slice of optimal strategy array
            optim_control_array_animallevel_popnpersp[:,vacc_strat_itr] = OptimThresholdVals_AnimalLevel_PopnPersp[:,percentile_slice_idx]
            optim_control_array_animallevel_indivpersp[:,vacc_strat_itr] = OptimThresholdVals_AnimalLevel_IndivPersp[:,percentile_slice_idx]

            #=
            Get the percentage of premises infected &
            vaccinated under the optimal strategy
            =#
            if compute_infections_vacc_flag == true

                # Based on scenario ID, pick out relevant batch number for the
                # optimal strategy within that scenario
                for rel_vacc_cost_itr = 1:n_rel_vacc_costs

                    # Animal level, population perspective
                    percentage_threshold_inf_exceeded_array_animallevel_popnpersp[rel_vacc_cost_itr,vacc_strat_itr,:],
                     percentage_infected_array_animallevel_popnpersp[rel_vacc_cost_itr,vacc_strat_itr,:],
                     percentage_vacc_array_animallevel_popnpersp[rel_vacc_cost_itr,vacc_strat_itr,:] =
                                    get_infection_vacc_outcomes(rel_vacc_cost_itr,
                                                                 vacc_strat_itr,
                                                                 batch_ID_offset,
                                                                 n_risk_measures,
                                                                 epi_outputs_agg_directory,
                                                                 optim_control_array_animallevel_popnpersp,
                                                                 prem_num,
                                                                 threshold_infections_vals,
                                                                 prctiles_infections_vaccs)

                    # Animal level, individual perspective
                    percentage_threshold_inf_exceeded_array_animallevel_indivpersp[rel_vacc_cost_itr,vacc_strat_itr,:],
                     percentage_infected_array_animallevel_indivpersp[rel_vacc_cost_itr,vacc_strat_itr,:],
                     percentage_vacc_array_animallevel_indivpersp[rel_vacc_cost_itr,vacc_strat_itr,:] =
                                    get_infection_vacc_outcomes(rel_vacc_cost_itr,
                                                                 vacc_strat_itr,
                                                                 batch_ID_offset,
                                                                 n_risk_measures,
                                                                 epi_outputs_agg_directory,
                                                                 optim_control_array_animallevel_indivpersp,
                                                                 prem_num,
                                                                 threshold_infections_vals,
                                                                 prctiles_infections_vaccs)

                end
            end
        end
        println("vacc_strat_itr: $vacc_strat_itr")
    end

    #--------------------------------------------------------------------------
    # Generate optimal strategy plots for each relative cost of vacc (in rel_cost_of_vacc_idxs)
    #--------------------------------------------------------------------------
    if plot_optim_strategy_flag == true
        for rel_vacc_cost_itr = 1:length(rel_cost_of_vacc_idxs)

            # Assign index corresponding to relative vaccine cost in this iteration
            # to variable
            val_cost_of_vacc_col_idx = rel_cost_of_vacc_idxs[rel_vacc_cost_itr]

            # Set title for the plot
            if val_cost_of_vacc_col_idx == 101
                rel_vacc_cost_val = 1
            elseif val_cost_of_vacc_col_idx ∈ collect(11:10:91)
                rel_vacc_cost_val = round((val_cost_of_vacc_col_idx*0.01)-0.01,digits=1)
            else
                rel_vacc_cost_val = round((val_cost_of_vacc_col_idx*0.01)-0.01,digits=2)
            end
            # title_string_popn_persp = "Population perspective, Cᵥ = $rel_vacc_cost_val"
            # title_string_indiv_persp = "Individual perspective, Cᵥ = $rel_vacc_cost_val"
            title_string_popn_persp = L"\textbf{Population \ perspective,} \, \textbf{C_V = %$rel_vacc_cost_val}"
            title_string_indiv_persp = L"\textbf{Individual \ perspective,} \, \textbf{C_V = %$rel_vacc_cost_val}"

            # Extract optimal strategy data
            optim_strat_vec_animallevel_popnpersp = optim_control_array_animallevel_popnpersp[val_cost_of_vacc_col_idx,:]
            optim_strat_vec_animallevel_indivpersp = optim_control_array_animallevel_indivpersp[val_cost_of_vacc_col_idx,:]

            #--  animallevel_popnpersp
            save_filename_point_plot_popn_persp = string("OptimBehaviour_SellkeRuns_TernaryPlots/animallevel_popn_persp_sellke_run_ternary_data_point_plot_", config_ID ,"_vacc_cost_ID_$(val_cost_of_vacc_col_idx)_julia.pdf")
            data_point_ternary(ternay_plot_pos_array,
                                    optim_strat_vec_animallevel_popnpersp,
                                    ternary_axes_labels,
                                    title_string_popn_persp,
                                    optim_strat_chosen_cmap,
                                    optim_strat_colourbar_label,
                                    optim_strat_cbar_range,
                                    optim_strat_colourbar_tick_pos,
                                    optim_strat_colourbar_tick_labels,
                                    label_fontsize,
                                    tick_label_fontsize,
                                    save_filename_point_plot_popn_persp,
                                    save_figure_flag)

            #--  animallevel_indpersp
            save_filename_point_plot_indiv_persp = string("OptimBehaviour_SellkeRuns_TernaryPlots/animallevel_indiv_persp_sellke_run_ternary_data_point_plot_", config_ID ,"_vacc_cost_ID_$(val_cost_of_vacc_col_idx)_julia.pdf")
            data_point_ternary(ternay_plot_pos_array,
                                    optim_strat_vec_animallevel_indivpersp,
                                    ternary_axes_labels,
                                    title_string_indiv_persp,
                                    optim_strat_chosen_cmap,
                                    optim_strat_colourbar_label,
                                    optim_strat_cbar_range,
                                    optim_strat_colourbar_tick_pos,
                                    optim_strat_colourbar_tick_labels,
                                    label_fontsize,
                                    tick_label_fontsize,
                                    save_filename_point_plot_indiv_persp,
                                    save_figure_flag)
        end
    end

    #--------------------------------------------------------------------------
    # Generate percentage of replicates exceeding a threshold amount of premises
    # infected under the optimal strategy with plots for each relative cost of vacc
    #--------------------------------------------------------------------------
    if plot_infections_vaccs_flag == true

       # Infection and vaccination stats not computed. Load from file instead
       if compute_infections_vacc_flag == false
            if config_ID == "cumbria"
                MAT_file = matopen("OptimBehaviour_InfVacc_PlotData/cumbria_inf_vacc_ternary_plot_data_julia.mat")
            elseif config_ID == "cumbria_alt"
                MAT_file = matopen("OptimBehaviour_InfVacc_PlotData/cumbria_alt_inf_vacc_ternary_plot_data_julia.mat")
            elseif config_ID == "devon"
                MAT_file = matopen("OptimBehaviour_InfVacc_PlotData/devon_vacc_ternary_plot_data_julia.mat")
            elseif config_ID == "devon_alt"
                MAT_file = matopen("OptimBehaviour_InfVacc_PlotData/devon_alt_vacc_ternary_plot_data_julia.mat")
            end
            percentage_threshold_inf_exceeded_array_animallevel_popnpersp = read(MAT_file, "percentage_threshold_inf_exceeded_array_animallevel_popnpersp")
            percentage_threshold_inf_exceeded_array_animallevel_indivpersp = read(MAT_file, "percentage_threshold_inf_exceeded_array_animallevel_indivpersp")
            close(MAT_file)
       end

       for rel_vacc_cost_itr = 1:length(rel_cost_of_vacc_idxs)

           # Assign index corresponding to relative vaccine cost in this iteration
           # to variable
           val_cost_of_vacc_col_idx = rel_cost_of_vacc_idxs[rel_vacc_cost_itr]

           # Set title for the plot
           if val_cost_of_vacc_col_idx == 101
               rel_vacc_cost_val = 1
           elseif val_cost_of_vacc_col_idx ∈ collect(11:10:91)
               rel_vacc_cost_val = round((val_cost_of_vacc_col_idx*0.01)-0.01,digits=1)
           else
               rel_vacc_cost_val = round((val_cost_of_vacc_col_idx*0.01)-0.01,digits=2)
           end
           # title_string_popn_persp = "Population perspective, Cᵥ = $rel_vacc_cost_val"
           # title_string_indiv_persp = "Individual perspective, Cᵥ = $rel_vacc_cost_val"
           title_string_popn_persp = L"\textbf{Population \ perspective,} \, \textbf{C_V = %$rel_vacc_cost_val}"
           title_string_indiv_persp = L"\textbf{Individual \ perspective,} \, \textbf{C_V = %$rel_vacc_cost_val}"

           # Extract percentage of replicates exceeding threshold under the optimal strategy data
           threshold_val_idx = 2 # To access threshold_infections_vals = [10 25 50 100 250 500 1000]
           exceed_threshold_infec_optim_strat_vec_animallevel_popnpersp = percentage_threshold_inf_exceeded_array_animallevel_popnpersp[val_cost_of_vacc_col_idx,:,threshold_val_idx]
           exceed_threshold_infec_optim_strat_vec_animallevel_indivpersp = percentage_threshold_inf_exceeded_array_animallevel_indivpersp[val_cost_of_vacc_col_idx,:,threshold_val_idx]

          #--  animallevel_popnpersp

           # Produce data point ternary plot
           save_filename_inf_point_plot_popn_persp = string("OptimBehaviour_InfecUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_popn_persp_sellke_run_exceed_infection_threshold_ternary_data_point_plot_", config_ID ,"_vacc_cost_ID_$(val_cost_of_vacc_col_idx)_julia.pdf")
           data_point_ternary(ternay_plot_pos_array,
                               exceed_threshold_infec_optim_strat_vec_animallevel_popnpersp,
                               ternary_axes_labels,
                               title_string_popn_persp ,
                               strat_exceed_threshold_infections_chosen_cmap,
                               strat_exceed_threshold_infections_colourbar_label,
                               strat_exceed_threshold_infections_cbar_range,
                               strat_exceed_threshold_infections_colourbar_tick_pos,
                               strat_exceed_threshold_infections_colourbar_tick_labels,
                               label_fontsize,
                               tick_label_fontsize,
                               save_filename_inf_point_plot_popn_persp,
                               save_figure_flag)

           #--  animallevel_indpersp

           # Produce data point ternary plot
           save_filename_inf_point_plot_indiv_persp = string("OptimBehaviour_InfecUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_indiv_persp_sellke_run_exceed_infection_threshold_ternary_data_point_plot_", config_ID ,"_vacc_cost_ID_$(val_cost_of_vacc_col_idx)_julia.pdf")
           data_point_ternary(ternay_plot_pos_array,
                                   exceed_threshold_infec_optim_strat_vec_animallevel_indivpersp,
                                   ternary_axes_labels,
                                   title_string_indiv_persp,
                                   strat_exceed_threshold_infections_chosen_cmap,
                                   strat_exceed_threshold_infections_colourbar_label,
                                   strat_exceed_threshold_infections_cbar_range,
                                   strat_exceed_threshold_infections_colourbar_tick_pos,
                                   strat_exceed_threshold_infections_colourbar_tick_labels,
                                   label_fontsize,
                                   tick_label_fontsize,
                                   save_filename_inf_point_plot_indiv_persp,
                                   save_figure_flag)
       end
    end

    #--------------------------------------------------------------------------
    # Generate vaccinations under the optimal strategy plots for each relative cost of vacc
    #--------------------------------------------------------------------------
    if plot_infections_vaccs_flag == true
        # Infection and vaccination stats not computed. Load from file instead
        if compute_infections_vacc_flag == false
             if config_ID == "cumbria"
                 MAT_file = matopen("OptimBehaviour_InfVacc_PlotData/cumbria_inf_vacc_ternary_plot_data_julia.mat")
             elseif config_ID == "cumbria_alt"
                 MAT_file = matopen("OptimBehaviour_InfVacc_PlotData/cumbria_alt_inf_vacc_ternary_plot_data_julia.mat")
             elseif config_ID == "devon"
                 MAT_file = matopen("OptimBehaviour_InfVacc_PlotData/devon_vacc_ternary_plot_data_julia.mat")
             elseif config_ID == "devon_alt"
                 MAT_file = matopen("OptimBehaviour_InfVacc_PlotData/devon_alt_vacc_ternary_plot_data_julia.mat")
             end
             percentage_vacc_array_animallevel_popnpersp = read(MAT_file, "percentage_vacc_array_animallevel_popnpersp")
             percentage_vacc_array_animallevel_indivpersp = read(MAT_file, "percentage_vacc_array_animallevel_indivpersp")
             close(MAT_file)
        end

        for rel_vacc_cost_itr = 1:length(rel_cost_of_vacc_idxs)

            # Assign index corresponding to relative vaccine cost in this iteration
            # to variable
            val_cost_of_vacc_col_idx = rel_cost_of_vacc_idxs[rel_vacc_cost_itr]

            # Set title for the plot
            if val_cost_of_vacc_col_idx == 101
                rel_vacc_cost_val = 1
            elseif val_cost_of_vacc_col_idx ∈ collect(11:10:91)
                rel_vacc_cost_val = round((val_cost_of_vacc_col_idx*0.01)-0.01,digits=1)
            else
                rel_vacc_cost_val = round((val_cost_of_vacc_col_idx*0.01)-0.01,digits=2)
            end
            # title_string_popn_persp = "Population perspective, Cᵥ = $rel_vacc_cost_val"
            # title_string_indiv_persp = "Individual perspective, Cᵥ = $rel_vacc_cost_val"
            title_string_popn_persp = L"\textbf{Population \ perspective,} \, \textbf{C_V = %$rel_vacc_cost_val}"
            title_string_indiv_persp = L"\textbf{Individual \ perspective,} \, \textbf{C_V = %$rel_vacc_cost_val}"

            # Extract median vaccinations under the optimal strategy data
            median_vaccs_optim_strat_vec_animallevel_popnpersp = percentage_vacc_array_animallevel_popnpersp[val_cost_of_vacc_col_idx,:,2]
            median_vaccs_optim_strat_vec_animallevel_indivpersp = percentage_vacc_array_animallevel_indivpersp[val_cost_of_vacc_col_idx,:,2]

            #--  animallevel_popnpersp

             # Produce data point ternary plot
             save_filename_vacc_point_plot_popn_persp = string("OptimBehaviour_VaccUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_popn_persp_sellke_run_median_vacc_ternary_data_point_plot_", config_ID ,"_vacc_cost_ID_$(val_cost_of_vacc_col_idx)_julia.pdf")
             data_point_ternary(ternay_plot_pos_array,
                                 median_vaccs_optim_strat_vec_animallevel_popnpersp,
                                 ternary_axes_labels,
                                 title_string_popn_persp,
                                 strat_vaccs_chosen_cmap,
                                 strat_vaccs_colourbar_label,
                                 strat_vaccs_cbar_range,
                                 strat_vaccs_colourbar_tick_pos,
                                 strat_vaccs_colourbar_tick_labels,
                                 label_fontsize,
                                 tick_label_fontsize,
                                 save_filename_vacc_point_plot_popn_persp,
                                 save_figure_flag)

             #--  animallevel_indpersp

             # Produce data point ternary plot
             save_filename_vacc_point_plot_indiv_persp = string("OptimBehaviour_VaccUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_indiv_persp_sellke_run_median_vacc_ternary_data_point_plot_", config_ID ,"_vacc_cost_ID_$(val_cost_of_vacc_col_idx)_julia.pdf")
             data_point_ternary(ternay_plot_pos_array,
                                     median_vaccs_optim_strat_vec_animallevel_indivpersp,
                                     ternary_axes_labels,
                                     title_string_popn_persp,
                                     strat_vaccs_chosen_cmap,
                                     strat_vaccs_colourbar_label,
                                     strat_vaccs_cbar_range,
                                     strat_vaccs_colourbar_tick_pos,
                                     strat_vaccs_colourbar_tick_labels,
                                     label_fontsize,
                                     tick_label_fontsize,
                                     save_filename_vacc_point_plot_indiv_persp,
                                     save_figure_flag)
        end
    end


    return percentage_threshold_inf_exceeded_array_animallevel_popnpersp::Array{Float64,3},
           percentage_threshold_inf_exceeded_array_animallevel_indivpersp::Array{Float64,3},
           percentage_infected_array_animallevel_popnpersp::Array{Float64,3},
           percentage_infected_array_animallevel_indivpersp::Array{Float64,3},
           percentage_vacc_array_animallevel_popnpersp::Array{Float64,3},
           percentage_vacc_array_animallevel_indivpersp::Array{Float64,3}
end

# SUPPORTING FUNCTIONS TO GENERATE TERNARY PLOT
"""
    data_point_ternary(ternay_plot_pos_array::Array{Float64,2},
                            data_vec::Array{Float64,1},
                            ternary_axes_labels::Array{String,1},
                            title_string::LaTeXString,
                            chosen_cmap::Any,
                            colourbar_label::String,
                            cbar_range::Tuple{Float64, Float64},
                            colourbar_tick_pos::Array{Float64,1},
                            colourbar_tick_labels::Array{String,1},
                            label_fontsize::Int64,
                            tick_label_fontsize::Int64,
                            save_filename_point_plot::String,
                            save_figure_flag::Bool)

Produce a ternary plot.

Inputs:
- `ternay_plot_pos_array::Array{Float64,2}`: The vaccine group compositions (three columns) for each scenario (rows).
- `data_vec::Array{Float64,1}`: Statistic to be plotted for each vaccine group scenario.
- `ternary_axes_labels::Array{String,1}`: Labels for each of the ternary axes.
- `title_string::LaTeXString`: Title for the ternary plot. Input "" for no title.
- `chosen_cmap::Any`: The colourmap applied to the markers.
- `colourbar_label::String`: Colourbar axis label.
- `cbar_range::Tuple{Float64, Float64}`: Range of the colourbar.
- `colourbar_tick_pos::Array{Float64,1}`: Set tick labels.
- `colourbar_tick_labels::Array{String,1}`: Set fontsize of colourbar tick labels
- `label_fontsize::Int64`: Set fontsize of axis labels.
- `tick_label_fontsize::Int64`: Set fontsize of tick labels.
- `save_filename_point_plot::String`: If saved, the filename plot is saved to.
- `save_figure_flag::Bool`: Flag to specify if plots should be saved to file

Outputs: None \n
Location: OptimThresholdVsCostRatio\\_SellkeRuns\\_GenerateTernaryPlots.jl
"""
function data_point_ternary(ternay_plot_pos_array::Array{Float64,2},
                            data_vec::Array{Float64,1},
                            ternary_axes_labels::Array{String,1},
                            title_string::LaTeXString,
                            chosen_cmap::Any,
                            colourbar_label::String,
                            cbar_range::Tuple{Float64, Float64},
                            colourbar_tick_pos::Array{Float64,1},
                            colourbar_tick_labels::Array{String,1},
                            label_fontsize::Int64,
                            tick_label_fontsize::Int64,
                            save_filename_point_plot::String,
                            save_figure_flag::Bool)

    # Initialise array to store ternary plot co-ords
    cartesian_coords = [zeros(eltype(ternay_plot_pos_array), size(ternay_plot_pos_array, 1)) zeros(eltype(ternay_plot_pos_array), size(ternay_plot_pos_array, 1))]

    # Populate ternary co-ords array
    for ii in 1:size(ternay_plot_pos_array,1)
        cartesian_coords[ii,:] = collect(tern2cart(ternay_plot_pos_array[ii,:]))'
    end

    # Generate ternary plot axis
    ternary_axes(title=title_string,
                   # titlefont = "Helvetica Bold",
                   xguide = ternary_axes_labels[2],
                   yguide = ternary_axes_labels[3],
                   zguide = ternary_axes_labels[1],
                   axis_arrows = false)
                   # xticks=[0,0.5,1],

    # Add markers
    p_ternary = scatter!(cartesian_coords[:,1],cartesian_coords[:,2],
                        marker_z = data_vec,
                        color = chosen_cmap,
                        clims = cbar_range,
                        colorbar_title = colourbar_label,
                        colorbar_titlefontsize = label_fontsize,
                        colorbar_ticks=[0, 1],
                        guidefont = 10,
                        markershape = :circle,
                        markersize = 10,
                        markeralpha = 1,
                        markerstrokewidth = 3,
                        markerstrokealpha = 1,
                        markerstrokecolor = :white,
                        markerstrokestyle = :solid,
                        legend=false)

    # Save the figure
    if save_figure_flag == true
        savefig(p_ternary, save_filename_point_plot)
    end

    return nothing
end

# SUPPORTING FUNCTIONS TO EXTRACT INFECTION & VACCINATION OUTCOMES
"""
    get_infection_vacc_outcomes(rel_vacc_cost_idx::Int64,
                                     vacc_strat_idx::Int64,
                                     batch_ID_offset::Int64,
                                     n_risk_measures::Int64,
                                     epi_outputs_agg_directory::String,
                                     optim_control_array::Array{Int64,2},
                                     prem_num::Int64,
                                     threshold_infections_vals::Array{Int64,1},
                                     prctiles_infections_vaccs::Array{Float64,1})

Computes infection and vaccination summary statistics to be used in ternary plots.

Inputs:
- `rel_cost_of_vacc_idxs::Array{Int64,1}`: Relative cost of vaccination to generate plots for
- `vacc_strat_idx::Int64`: ID for the vaccine group composition scenario
- `batch_ID_offset::Int64`: Added to the batch_ID (for accessing model output files)
- `n_risk_measures::Int64`: Number of different reactive control strats
- `epi_outputs_agg_directory::String`: Directory where model outputs reside
- `optim_control_array::Array{Int64,2}`: Optimal strategy for each combination of relative vaccine cost (rows) and vaccine group scenario (columns)
- `prem_num::Int64`: Number of premises in the landscape
- `threshold_infections_vals::Array{Int64,1}`: pexceedence of threshold infection to be computed
- `prctiles_infections_vaccs::Array{Float64,1}`: percentile values to be computed

Outputs: Percentages for replicates exceeding infection threshold, premises infected & vaccinated under the optimal strategy
- `percentage_threshold_inf_exceeded_vec::Array{Float64,1}`
- `percentage_infected_vec::Array{Float64,1}`
- `percentage_vacc_vec::Array{Float64,1}`

Location: OptimThresholdVsCostRatio\\_SellkeRuns\\_GenerateTernaryPlots.jl
"""
function get_infection_vacc_outcomes(rel_vacc_cost_idx::Int64,
                                     vacc_strat_idx::Int64,
                                     batch_ID_offset::Int64,
                                     n_risk_measures::Int64,
                                     epi_outputs_agg_directory::String,
                                     optim_control_array::Array{Int64,2},
                                     prem_num::Int64,
                                     threshold_infections_vals::Array{Int64,1},
                                     prctiles_infections_vaccs::Array{Float64,1})

    # Premises level, population perspective
    optim_stratID = optim_control_array[rel_vacc_cost_idx,vacc_strat_idx] + 1
        # optim_control_array(rel_vacc_cost_itr,vacc_strat_itr)
        # returns an integer value [0,1,,n_risk_measures].
        # Add 1 to get the configuration within the scenario
        # set to be used.

    # Chosen batch ID
    chosen_batch_ID = batch_ID_offset + ((vacc_strat_idx-1)*n_risk_measures) + optim_stratID

    # Load the data from the relevant Batch ID
    optim_strat_epi_outputs_filename = string(epi_outputs_agg_directory,"PremPerDiseaseState_BatchID$(chosen_batch_ID).txt")
    optim_strat_epi_outputs_array = readdlm(optim_strat_epi_outputs_filename)

    # For that batch ID, get the number of premises infected and vaccinated
    optim_strat_infections = optim_strat_epi_outputs_array[:,6]
    optim_strat_vaccs = optim_strat_epi_outputs_array[:,7]

    # Compute the number of replicates exceeding the specified infected
    # premises threshold values
    percentage_threshold_inf_exceeded_vec = zeros(Float64,length(threshold_infections_vals))
    for threshold_val_itr = 1:length(threshold_infections_vals)
        percentage_threshold_inf_exceeded_vec[threshold_val_itr] = (sum(optim_strat_infections .>= threshold_infections_vals[threshold_val_itr])/length(optim_strat_infections))*100
    end

    # Compute requested percentiles for infection data across all replicates
    # Turn into percentage based on total number of premises
    n_prctiles_infections_vaccs = length(prctiles_infections_vaccs)
    percentage_infected_vec = zeros(Float64,n_prctiles_infections_vaccs)
    percentage_vacc_vec = zeros(Float64,n_prctiles_infections_vaccs)
    for prctile_itr = 1:n_prctiles_infections_vaccs
        percentage_infected_vec[prctile_itr] = (percentile(optim_strat_infections,prctiles_infections_vaccs[prctile_itr])/prem_num)*100
        percentage_vacc_vec[prctile_itr] = (percentile(optim_strat_vaccs,prctiles_infections_vaccs[prctile_itr])/prem_num)*100
    end

    return percentage_threshold_inf_exceeded_vec::Array{Float64,1},
                percentage_infected_vec::Array{Float64,1},
                percentage_vacc_vec::Array{Float64,1}
end
