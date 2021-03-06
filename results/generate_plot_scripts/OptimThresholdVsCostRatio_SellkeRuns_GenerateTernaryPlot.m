%Purpose:
% Produce outputs of optimal strategy (out of a collection of strategies)
% based on three different groups
% (i) Initial control adopters (enact control from outset of simualation)
% (ii) Reactive responders (implement control based on epidemic situation)
% (iii) Never implement control
%
% Plot produced for a single relative cost of vaccination
% Outbreak simulations used the Sellke construction.
%
% Also can produce ternary plots for:
% - cost under optimal strategy
% - estimated infections (for specified percentiles) under optimal strategy
% - estimated vaccinations (for specified percentiles) under optimal strategy
%
% Cite Ternary Plots package as
% Ulrich Theune (2021). Ternary Plots (https://www.mathworks.com/matlabcentral/fileexchange/7210-ternary-plots),
% MATLAB Central File Exchange. Retrieved October 4, 2021.
%
% MATLAB version: R2021b
% Date: 3rd November 2021
%--------------------------------------------------------------------------


function [percentage_threshold_inf_exceeded_array_animallevel_popnpersp,...
          percentage_threshold_inf_exceeded_array_animallevel_indivpersp,...
          percentage_infected_array_animallevel_popnpersp,...
          percentage_infected_array_animallevel_indivpersp,...
          percentage_vacc_array_animallevel_popnpersp,...
          percentage_vacc_array_animallevel_indivpersp] = OptimThresholdVsCostRatio_SellkeRuns_GenerateTernaryPlot(save_figure_flag,...
                                                                    config_ID,...
                                                                    scen_IDs,...
                                                                    batch_ID_offset ,...
                                                                    rel_cost_of_vacc_idxs,...
                                                                    compute_optim_strategy_flag,...
                                                                    plot_optim_strategy_flag,...
                                                                    plot_costs_flag,...
                                                                    plot_infections_vaccs_flag,...
                                                                    compute_infections_vacc_flag,...
                                                                    prctiles_infections_vaccs,...
                                                                    threshold_infections_vals)

    %%
    %--------------------------------------------------------------------------
    %%% Set ternary plot variables
    %--------------------------------------------------------------------------

    % Specify number of control strategies considered
    n_risk_measures = 11; % 0km, 1km, .... , 9km, 10km distance based measures

    % Axes labels
    ternary_axes_labels = {'Precautionary (%)','Reactionary (%)','Non-vaccinators (%)'};

    % Font sizes
    label_fontsize = 22;
    tick_label_fontsize = 20;

    % Colourbar
    cbar_range = [-0.5 n_risk_measures-0.5];
    colourbar_tick_pos = 0:1:(n_risk_measures-1);
    colourbar_tick_labels = {'0','1','2','3','4','5','6','7','8','9','10'};
    colourbar_label = 'Risk threshold: Notified premises distance (km)';

    % Setup colour map options for optimal strategy plots
    CT_optim_strat = cbrewer('seq', 'Reds', n_risk_measures);

    my_colourmap = [0.5 0.5 0.5;
                    CT_optim_strat(end:-1:2,:)]; % Have first entry as grey

    chosen_cmap = my_colourmap;

    % Setup colourbar options for cost under optimal strategy plots
    CT_blue_cost_plots = cbrewer('seq', 'Blues', 1000);
    strat_cost_chosen_cmap = CT_blue_cost_plots;
    strat_cost_colourbar_label_popn_persp = 'Population perspective cost (using optimal strategy)';
    strat_cost_colourbar_label_indiv_persp = 'Individual perspective cost (using optimal strategy)';
    strat_cost_cbar_range = [0 1500000];
    strat_cost_colourbar_tick_pos = [0 250000 500000 750000 1000000 1250000 1500000];
    strat_cost_colourbar_tick_labels = {'0','250000','500000','750000','1000000','1250000','1500000'};

    % Setup colourbar options for exceeding threshold infection value under optimal strategy plots
    strat_exceed_threshold_infections_chosen_cmap = cbrewer('seq', 'Oranges', 1000);
    strat_exceed_threshold_infections_colourbar_label = 'Replicates with 25+ premises infected (%)';
    strat_exceed_threshold_infections_cbar_range = [0 100];
    strat_exceed_threshold_infections_colourbar_tick_pos = 0:10:100;
    strat_exceed_threshold_infections_colourbar_tick_labels = 0:10:100;

    % Setup colourbar options for vaccinations under optimal strategy plots
    strat_vaccs_chosen_cmap = cbrewer('seq', 'Blues', 1000);
    strat_vaccs_colourbar_label = 'Premises vaccinated (%)';
    strat_vaccs_cbar_range = [0 100];
    strat_vaccs_colourbar_tick_pos = [0 10 20 30 40 50 60 70 80 90 100];
    strat_vaccs_colourbar_tick_labels = [0 10 20 30 40 50 60 70 80 90 100];

    % Specify amount of infection related thresholds in use
    n_prctiles_infections_vaccs = length(prctiles_infections_vaccs);
    n_threshold_infections_vals = length(threshold_infections_vals);

    %%
    %--------------------------------------------------------------------------
    %%% Load optimum strategy data
    %--------------------------------------------------------------------------

    % Get number of scenarios
    n_scens = length(scen_IDs);

    % Set aggregated epidemiological output file directory to be accessed
    % and number of premises in use
    if strcmp(config_ID,'cumbria') || strcmp(config_ID,'cumbria_alt')
        epi_outputs_agg_directory = '../GBCountyModelSimnOutputs/Cumbria_EpiOutputs_Aggregated/';
        prem_num = 3784;
    elseif strcmp(config_ID,'devon') || strcmp(config_ID,'devon_alt')
        epi_outputs_agg_directory = '../GBCountyModelSimnOutputs/Devon_EpiOutputs_Aggregated/';
        prem_num = 5343;
    else
        error('Invalid config_ID provided.')
    end

    % Initialise filename vector
    filename_suffix_vec = cell(n_scens,1);

    % Populate filename vector
    for scen_itr = 1:n_scens

        % Get scenario ID for current iteration
        current_itr_scen_ID = scen_IDs(scen_itr);

        % Assign filename to cel
        if strcmp(config_ID,'cumbria')
            filename_suffix_vec{scen_itr} = ['Cumbria_vacc_distance_risk_measure_scenID', num2str(current_itr_scen_ID)];
        elseif strcmp(config_ID,'cumbria_alt')
            filename_suffix_vec{scen_itr} = ['Cumbria_alt_pathogen_vacc_distance_risk_measure_scenID', num2str(current_itr_scen_ID)];
        elseif strcmp(config_ID,'devon')
            filename_suffix_vec{scen_itr} = ['Devon_vacc_distance_risk_measure_scenID', num2str(current_itr_scen_ID)];
        elseif strcmp(config_ID,'devon_alt')
            filename_suffix_vec{scen_itr} = ['Devon_alt_pathogen_vacc_distance_risk_measure_scenID', num2str(current_itr_scen_ID)];
        else
            error('Invalid config_ID provided.')
        end
    end

    % Assign number of vaccine group combinations tested to variable
    n_vacc_uptake_strats = length(filename_suffix_vec);

    %%
    %--------------------------------------------------------------------------
    %%% Construct ternary plot position data
    %--------------------------------------------------------------------------

    % Position array format: Column 1: Initial adopters fraction
    %                        Column 2: Reactioners fraction
    %                        Column 3: Non-vaccinators fraction
    ternay_plot_pos_array = zeros(n_scens,3);
    % ternay_plot_pos_array = [0 1 0;
    %                          0.05 0.95 0;
    %                          0.1 0.9 0;
    %                          0.25 0.75 0;
    %                          0.5 0.5 0];

    % Load group composition per scenario from data file
    control_param_input_args_array = readmatrix('../../src/spatial_simns/run_model_scripts/parameter_combination_files/generic_model_control_param_val_array.txt');
        % Three column array.
        % Column 1: Proportion in "Initial adopters" group
        % Column 2: Proportion inf "Reactioners" group
        % Column 3: Risk threshold measure (NOT USED HERE)

    % Assign "Initial adopters" & "Reactioners" group columns to array
    % Take every "n_risk_measures"th row from control_param_input_args_array,
    % as that corresponds to moving to the next scenario chunk that used a
    % different group composition
    ternay_plot_pos_array(:,1:2) = control_param_input_args_array(1:n_risk_measures:end,1:2);

    % Populate "No control" group using data on other groups
    ternay_plot_pos_array(:,3) = ones(n_scens,1) - sum(ternay_plot_pos_array(:,1:2),2);

    % Error check
    if sum((sum(ternay_plot_pos_array,2) ~= 1)) > 0
        error('Rows in ternay_plot_pos_array do not sum to 1 as expected.')
    end

    %%
    %--------------------------------------------------------------------------
    %%% Construct optimal strategy arrays for all relative vaccine costs
    %%% Assign the cost under the optimal strategy to array
    %--------------------------------------------------------------------------

    % Specify the number of relative vaccination costs tested
    n_rel_vacc_costs = 101;

    % Specify the percentile to be used in plots
    % Slice idx: 1 - mean; 2 - 2.5th percentile; 3 - median; 4 - 97.5th
    % percentile; 5 - minimum; 6 - maximum.
    percentile_slice_idx = 3;

    % Initialise storage arrays for the optimal strategy
    optim_control_array_premlevel_popnpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats);
    optim_control_array_premlevel_indivpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats);

    optim_control_array_animallevel_popnpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats);
    optim_control_array_animallevel_indivpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats);

    % Initialise storage arrays to record the cost under the strategy
    % identified as optimal
    cost_of_optim_control_array_premlevel_popnpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats);
    cost_of_optim_control_array_premlevel_indivpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats);

    cost_of_optim_control_array_animallevel_popnpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats);
    cost_of_optim_control_array_animallevel_indivpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats);

    % Initialise storage arrays to record the percentage of replicates that
    % exceeded a specified number of premises infected under the strategy identified as optimal
    percentage_threshold_inf_exceeded_array_animallevel_popnpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats,n_threshold_infections_vals);
    percentage_threshold_inf_exceeded_array_animallevel_indivpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats,n_threshold_infections_vals);

    % Initialise storage arrays to record the percentage of premises infected
    % under the strategy identified as optimal
    percentage_infected_array_animallevel_popnpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats,n_prctiles_infections_vaccs);
    percentage_infected_array_animallevel_indivpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats,n_prctiles_infections_vaccs);

    % Initialise storage arrays to record the percentage of premises vaccinated
    % under the strategy identified as optimal
    percentage_vacc_array_animallevel_popnpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats,n_prctiles_infections_vaccs);
    percentage_vacc_array_animallevel_indivpersp = zeros(n_rel_vacc_costs,n_vacc_uptake_strats,n_prctiles_infections_vaccs);

    % Populate arrays, iterate over each scenario
    for vacc_strat_itr = 1:n_vacc_uptake_strats

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Get the optimal strategy data %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if compute_optim_strategy_flag == true
            % Set relevant data file names
            OptimThresholdVals_PremLevel_PopnPersp_FileName = ['OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_PremLevel_PopnPersp_' filename_suffix_vec{vacc_strat_itr} '.txt'];
            OptimThresholdVals_AnimalLevel_PopnPersp_FileName = ['OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_PopnPersp_' filename_suffix_vec{vacc_strat_itr} '.txt'];

            OptimThresholdVals_PremLevel_IndivPersp_FileName = ['OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_PremLevel_IndivPersp_' filename_suffix_vec{vacc_strat_itr} '.txt']';
            OptimThresholdVals_AnimalLevel_IndivPersp_FileName = ['OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_IndivPersp_' filename_suffix_vec{vacc_strat_itr} '.txt'];

            % Load optimal strategy data
            OptimThresholdVals_PremLevel_PopnPersp = readmatrix(OptimThresholdVals_PremLevel_PopnPersp_FileName);
            OptimThresholdVals_AnimalLevel_PopnPersp = readmatrix(OptimThresholdVals_AnimalLevel_PopnPersp_FileName);

            OptimThresholdVals_PremLevel_IndivPersp = readmatrix(OptimThresholdVals_PremLevel_IndivPersp_FileName);
            OptimThresholdVals_AnimalLevel_IndivPersp = readmatrix(OptimThresholdVals_AnimalLevel_IndivPersp_FileName);

            % Pick out relevant slice of optimal strategy array
            optim_control_array_premlevel_popnpersp(:,vacc_strat_itr) = OptimThresholdVals_PremLevel_PopnPersp(:,percentile_slice_idx);
            optim_control_array_premlevel_indivpersp(:,vacc_strat_itr) = OptimThresholdVals_PremLevel_IndivPersp(:,percentile_slice_idx);
            optim_control_array_animallevel_popnpersp(:,vacc_strat_itr) = OptimThresholdVals_AnimalLevel_PopnPersp(:,percentile_slice_idx);
            optim_control_array_animallevel_indivpersp(:,vacc_strat_itr) = OptimThresholdVals_AnimalLevel_IndivPersp(:,percentile_slice_idx);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Get the cost under the optimal strategy %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_costs_flag == true

                % Set relevant data file names
                optim_strat_cost_filename = ['OptimBehaviour_SummStatsTotalCostsOutputs/OptimBehav_SummStatTotalCostData_' filename_suffix_vec{vacc_strat_itr} '.mat'];
                optim_strat_cost_data_array = load(optim_strat_cost_filename);

                % Load optimal strategy data
                OptimStratCost_PremLevel_PopnPersp = optim_strat_cost_data_array.SummStatTotalCost_PremLevel_PopnPersp;
                OptimStratCost_AnimalLevel_PopnPersp = optim_strat_cost_data_array.SummStatTotalCost_AnimalLevel_PopnPersp;

                OptimStratCost_PremLevel_IndivPersp = optim_strat_cost_data_array.SummStatTotalCost_PremLevel_IndivPersp;
                OptimStratCost_AnimalLevel_IndivPersp = optim_strat_cost_data_array.SummStatTotalCost_AnimalLevel_IndivPersp;

                % Pick out cost associated with the strategy identified as optimal for
                % a given relative cost of vacciantion
                for rel_vacc_cost_itr = 1:n_rel_vacc_costs

                    % Premises level, population perspective
                    optim_stratID_premlevel_popnpersp = optim_control_array_premlevel_popnpersp(rel_vacc_cost_itr,vacc_strat_itr) + 1;
                    cost_of_optim_control_array_premlevel_popnpersp(rel_vacc_cost_itr,vacc_strat_itr) = OptimStratCost_PremLevel_PopnPersp(rel_vacc_cost_itr,optim_stratID_premlevel_popnpersp,percentile_slice_idx);
                        % optim_control_array_premlevel_popnpersp(rel_vacc_cost_itr,vacc_strat_itr)
                        % returns an integer value [0,1,...,9,10].
                        % Corresponds with columns [1,2,...,10,11] of OptimStratCost_PremLevel_PopnPersp,
                        % hence the +1

                    % Premises level, individual perspective
                    optim_stratID_premlevel_indivpersp = optim_control_array_premlevel_indivpersp(rel_vacc_cost_itr,vacc_strat_itr) + 1;
                    cost_of_optim_control_array_premlevel_indivpersp(rel_vacc_cost_itr,vacc_strat_itr) = OptimStratCost_PremLevel_IndivPersp(rel_vacc_cost_itr,optim_stratID_premlevel_indivpersp,percentile_slice_idx);

                    % Animal level, population perspective
                    optim_stratID_animallevel_popnpersp = optim_control_array_animallevel_popnpersp(rel_vacc_cost_itr,vacc_strat_itr) + 1;
                    cost_of_optim_control_array_animallevel_popnpersp(rel_vacc_cost_itr,vacc_strat_itr) = OptimStratCost_AnimalLevel_PopnPersp(rel_vacc_cost_itr,optim_stratID_animallevel_popnpersp,percentile_slice_idx);

                    % Animal level, individual perspective
                    optim_stratID_animallevel_indivpersp = optim_control_array_animallevel_indivpersp(rel_vacc_cost_itr,vacc_strat_itr) + 1;
                    cost_of_optim_control_array_animallevel_indivpersp(rel_vacc_cost_itr,vacc_strat_itr) = OptimStratCost_AnimalLevel_IndivPersp(rel_vacc_cost_itr,optim_stratID_animallevel_indivpersp,percentile_slice_idx);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Get the percentage of premises infected & %%%%%%%%%
            %%% vaccinated under the optimal strategy    %%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if compute_infections_vacc_flag == true

                % Based on scenario ID, pick out relevant batch number for the
                % optimal strategy within that scenario
                for rel_vacc_cost_itr = 1:n_rel_vacc_costs

                    % Animal level, population perspective
                    [percentage_threshold_inf_exceeded_array_animallevel_popnpersp(rel_vacc_cost_itr,vacc_strat_itr,:),...
                     percentage_infected_array_animallevel_popnpersp(rel_vacc_cost_itr,vacc_strat_itr,:),...
                     percentage_vacc_array_animallevel_popnpersp(rel_vacc_cost_itr,vacc_strat_itr,:)] =...
                                    get_infection_vacc_outcomes(rel_vacc_cost_itr,...
                                                                 vacc_strat_itr,...
                                                                 batch_ID_offset,...
                                                                 n_risk_measures,...
                                                                 epi_outputs_agg_directory,...
                                                                 optim_control_array_animallevel_popnpersp,...
                                                                 prem_num,...
                                                                 threshold_infections_vals,...
                                                                 prctiles_infections_vaccs);

                    % Animal level, individual perspective
                    [percentage_threshold_inf_exceeded_array_animallevel_indivpersp(rel_vacc_cost_itr,vacc_strat_itr,:),...
                     percentage_infected_array_animallevel_indivpersp(rel_vacc_cost_itr,vacc_strat_itr,:),...
                     percentage_vacc_array_animallevel_indivpersp(rel_vacc_cost_itr,vacc_strat_itr,:)] =...
                                    get_infection_vacc_outcomes(rel_vacc_cost_itr,...
                                                                 vacc_strat_itr,...
                                                                 batch_ID_offset,...
                                                                 n_risk_measures,...
                                                                 epi_outputs_agg_directory,...
                                                                 optim_control_array_animallevel_indivpersp,...
                                                                 prem_num,...
                                                                 threshold_infections_vals,...
                                                                 prctiles_infections_vaccs);
                end
            end
        end
        fprintf('vacc_strat_itr: %d \n',vacc_strat_itr)
    end

    %%
    %--------------------------------------------------------------------------
    %%% Map to plot values (optimal strategy arrays for all relative vaccine costs)
    %--------------------------------------------------------------------------

    % Map strategy distance threshold original value to colourbar value
    original_vals = 0:1:10;
    plot_vals = 0:1:10;
    % original_vals = [0,0.5,1,3,5,10];
    % plot_vals = [0,1,2,3,4,5];
    for alter_val_idx = 1:length(original_vals)
        current_itr_original_val = original_vals(alter_val_idx);
        current_itr_plot_val = plot_vals(alter_val_idx);

        optim_control_array_premlevel_popnpersp = modify_value(optim_control_array_premlevel_popnpersp,current_itr_original_val,current_itr_plot_val);
        optim_control_array_premlevel_indivpersp = modify_value(optim_control_array_premlevel_indivpersp,current_itr_original_val,current_itr_plot_val);

        optim_control_array_animallevel_popnpersp = modify_value(optim_control_array_animallevel_popnpersp,current_itr_original_val,current_itr_plot_val);
        optim_control_array_animallevel_indivpersp = modify_value(optim_control_array_animallevel_indivpersp,current_itr_original_val,current_itr_plot_val);
    end

    %%
    %--------------------------------------------------------------------------
    %%% Generate optimal strategy plots for each relative cost of vacc (in rel_cost_of_vacc_idxs)
    %--------------------------------------------------------------------------
    if plot_optim_strategy_flag == true
        optim_strat_ternary_flag = true;
        for rel_vacc_cost_itr = 1:length(rel_cost_of_vacc_idxs)

            % Assign index corresponding to relative vaccine cost in this iteration
            % to variable
            val_cost_of_vacc_col_idx = rel_cost_of_vacc_idxs(rel_vacc_cost_itr);

            % Set title for the plot
            rel_vacc_cost_val = (val_cost_of_vacc_col_idx*0.01)-0.01;
            title_string_popn_persp = ['Population perspective, C_{V} = ', num2str(rel_vacc_cost_val)];
            title_string_indiv_persp = ['Individual perspective, C_{V} = ', num2str(rel_vacc_cost_val)];

            % Extract optimal strategy data
            optim_strat_vec_animallevel_popnpersp = optim_control_array_animallevel_popnpersp(val_cost_of_vacc_col_idx,:);
            optim_strat_vec_animallevel_indivpersp = optim_control_array_animallevel_indivpersp(val_cost_of_vacc_col_idx,:);


            %--  animallevel_popnpersp

            % Produce data point ternary plot
            savefile_name_point_plot = ['OptimBehaviour_SellkeRuns_TernaryPlots/animallevel_popn_persp_sellke_run_ternary_data_point_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx),'.pdf'];
            data_point_ternary(ternay_plot_pos_array,...
                                    optim_strat_vec_animallevel_popnpersp,...
                                    optim_strat_ternary_flag,...
                                    ternary_axes_labels,...
                                    title_string_popn_persp,...
                                    chosen_cmap,...
                                    colourbar_label,...
                                    cbar_range,...
                                    colourbar_tick_pos,...
                                    colourbar_tick_labels,...
                                    label_fontsize,...
                                    tick_label_fontsize,...
                                    savefile_name_point_plot,...
                                    save_figure_flag)

            %--  animallevel_indpersp

            % Produce data point ternary plot
            savefile_name_point_plot = ['OptimBehaviour_SellkeRuns_TernaryPlots/animallevel_indiv_persp_sellke_run_ternary_data_point_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx),'.pdf'];
            data_point_ternary(ternay_plot_pos_array,...
                                    optim_strat_vec_animallevel_indivpersp,...
                                    optim_strat_ternary_flag,...
                                    ternary_axes_labels,...
                                    title_string_indiv_persp,...
                                    chosen_cmap,...
                                    colourbar_label,...
                                    cbar_range,...
                                    colourbar_tick_pos,...
                                    colourbar_tick_labels,...
                                    label_fontsize,...
                                    tick_label_fontsize,...
                                    savefile_name_point_plot,...
                                    save_figure_flag)
        end
    end
    %%
    %--------------------------------------------------------------------------
    %%% Generate the cost under the optimal strategy plots for each relative cost of vacc
    %--------------------------------------------------------------------------
    if plot_costs_flag == true
        optim_strat_ternary_flag = false;
        for rel_vacc_cost_itr = 1:length(rel_cost_of_vacc_idxs)

            % Assign index corresponding to relative vaccine cost in this iteration
            % to variable
            val_cost_of_vacc_col_idx = rel_cost_of_vacc_idxs(rel_vacc_cost_itr);

            % Set title for the plot
            rel_vacc_cost_val = (val_cost_of_vacc_col_idx*0.01)-0.01;
            % title_string = ['Relative cost of vaccination: ', num2str(rel_vacc_cost_val)];
            title_string_popn_persp = ['Population perspective, C_{V} = ', num2str(rel_vacc_cost_val)];
            title_string_indiv_persp = ['Individual perspective, C_{V} = ', num2str(rel_vacc_cost_val)];

            % Extract cost under the optimal strategy data
%             cost_of_optim_strat_vec_premlevel_popnpersp = cost_of_optim_control_array_premlevel_popnpersp(val_cost_of_vacc_col_idx,:);
%             cost_of_optim_strat_vec_premlevel_indivpersp = cost_of_optim_control_array_premlevel_indivpersp(val_cost_of_vacc_col_idx,:);
            cost_of_optim_strat_vec_animallevel_popnpersp = cost_of_optim_control_array_animallevel_popnpersp(val_cost_of_vacc_col_idx,:);
            cost_of_optim_strat_vec_animallevel_indivpersp = cost_of_optim_control_array_animallevel_indivpersp(val_cost_of_vacc_col_idx,:);

            %--  animallevel_popnpersp

             % Produce surface colour ternary plot
            savefile_name_surface_plot = ['OptimBehaviour_CostUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_popn_persp_sellke_run_strategy_cost_ternary_surface_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx)];
            filled_colour_ternary(ternay_plot_pos_array,...
                                        cost_of_optim_strat_vec_animallevel_popnpersp,...
                                        optim_strat_ternary_flag,...
                                        ternary_axes_labels,...
                                        title_string_popn_persp,...
                                        strat_cost_chosen_cmap,...
                                        strat_cost_colourbar_label_popn_persp,...
                                        strat_cost_cbar_range,...
                                        strat_cost_colourbar_tick_pos,...
                                        strat_cost_colourbar_tick_labels,...
                                        label_fontsize,...
                                        tick_label_fontsize,...
                                        savefile_name_surface_plot,...
                                        save_figure_flag)

            % Produce data point ternary plot
            savefile_name_point_plot = ['OptimBehaviour_CostUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_popn_persp_sellke_run_strategy_cost_ternary_data_point_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx)];
            data_point_ternary(ternay_plot_pos_array,...
                                    cost_of_optim_strat_vec_animallevel_popnpersp,...
                                    optim_strat_ternary_flag,...
                                    ternary_axes_labels,...
                                    title_string_popn_persp,...
                                    strat_cost_chosen_cmap,...
                                    strat_cost_colourbar_label_popn_persp,...
                                    strat_cost_cbar_range,...
                                    strat_cost_colourbar_tick_pos,...
                                    strat_cost_colourbar_tick_labels,...
                                    label_fontsize,...
                                    tick_label_fontsize,...
                                    savefile_name_point_plot,...
                                    save_figure_flag)

            %--  animallevel_indpersp

            % Produce data point ternary plot
            savefile_name_point_plot = ['OptimBehaviour_CostUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_indiv_persp_sellke_run_strategy_cost_ternary_data_point_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx)];
            data_point_ternary(ternay_plot_pos_array,...
                                    cost_of_optim_strat_vec_animallevel_indivpersp,...
                                    ternary_axes_labels,...
                                    title_string_indiv_persp,...
                                    strat_cost_chosen_cmap,...
                                    strat_cost_colourbar_label_indiv_persp,...
                                    strat_cost_cbar_range,...
                                    strat_cost_colourbar_tick_pos,...
                                    strat_cost_colourbar_tick_labels,...
                                    label_fontsize,...
                                    tick_label_fontsize,...
                                    savefile_name_point_plot,...
                                    save_figure_flag)
        end
    end

    %%
    %--------------------------------------------------------------------------
    %%% Generate percentage of replicates exceeding a threshold amount of premises
    %%% infected under the optimal strategy with plots for each relative cost of vacc
    %--------------------------------------------------------------------------
    if plot_infections_vaccs_flag == true
        optim_strat_ternary_flag = false;

        % Infection and vaccination stats not computed. Load from file instead
        if compute_infections_vacc_flag == false
            if strcmp(config_ID,'cumbria')
                load('OptimBehaviour_InfVacc_PlotData/cumbria_inf_vacc_ternary_plot_data.mat',...
                    'percentage_threshold_inf_exceeded_array_animallevel_popnpersp',...
                    'percentage_threshold_inf_exceeded_array_animallevel_indivpersp')
            elseif strcmp(config_ID,'cumbria_alt')
                load('OptimBehaviour_InfVacc_PlotData/cumbria_alt_inf_vacc_ternary_plot_data.mat',...
                    'percentage_threshold_inf_exceeded_array_animallevel_popnpersp',...
                    'percentage_threshold_inf_exceeded_array_animallevel_indivpersp')
            elseif strcmp(config_ID,'devon')
                load('OptimBehaviour_InfVacc_PlotData/devon_inf_vacc_ternary_plot_data.mat',...
                    'percentage_threshold_inf_exceeded_array_animallevel_popnpersp',...
                    'percentage_threshold_inf_exceeded_array_animallevel_indivpersp')
            elseif strcmp(config_ID,'devon_alt')
                load('OptimBehaviour_InfVacc_PlotData/devon_alt_inf_vacc_ternary_plot_data.mat',...
                    'percentage_threshold_inf_exceeded_array_animallevel_popnpersp',...
                    'percentage_threshold_inf_exceeded_array_animallevel_indivpersp')
            end
        end

        for rel_vacc_cost_itr = 1:length(rel_cost_of_vacc_idxs)

            % Assign index corresponding to relative vaccine cost in this iteration
            % to variable
            val_cost_of_vacc_col_idx = rel_cost_of_vacc_idxs(rel_vacc_cost_itr);

            % Set title for the plot
            rel_vacc_cost_val = (val_cost_of_vacc_col_idx*0.01)-0.01;
            title_string_popn_persp = ['Population perspective, C_{V} = ', num2str(rel_vacc_cost_val)];
            title_string_indiv_persp = ['Individual perspective, C_{V} = ', num2str(rel_vacc_cost_val)];

            % Extract percentage of replicates exceeding threshold under the optimal strategy data
            threshold_val_idx = 2; % To access threshold_infections_vals = [10 25 50 100 250 500 1000]
            exceed_threshold_infec_optim_strat_vec_animallevel_popnpersp = percentage_threshold_inf_exceeded_array_animallevel_popnpersp(val_cost_of_vacc_col_idx,:,threshold_val_idx);
            exceed_threshold_infec_optim_strat_vec_animallevel_indivpersp = percentage_threshold_inf_exceeded_array_animallevel_indivpersp(val_cost_of_vacc_col_idx,:,threshold_val_idx);

            %--  animallevel_popnpersp

            % Produce data point ternary plot
            savefile_name_point_plot = ['OptimBehaviour_InfecUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_popn_persp_sellke_run_exceed_infection_threshold_ternary_data_point_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx),'.pdf'];
            data_point_ternary(ternay_plot_pos_array,...
                                    exceed_threshold_infec_optim_strat_vec_animallevel_popnpersp,...
                                    optim_strat_ternary_flag,...
                                    ternary_axes_labels,...
                                    title_string_popn_persp ,...
                                    strat_exceed_threshold_infections_chosen_cmap,...
                                    strat_exceed_threshold_infections_colourbar_label,...
                                    strat_exceed_threshold_infections_cbar_range,...
                                    strat_exceed_threshold_infections_colourbar_tick_pos,...
                                    strat_exceed_threshold_infections_colourbar_tick_labels,...
                                    label_fontsize,...
                                    tick_label_fontsize,...
                                    savefile_name_point_plot,...
                                    save_figure_flag)

            %--  animallevel_indpersp

            % Produce data point ternary plot
            savefile_name_point_plot = ['OptimBehaviour_InfecUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_indiv_persp_sellke_run_exceed_infection_threshold_ternary_data_point_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx),'.pdf'];
            data_point_ternary(ternay_plot_pos_array,...
                                    exceed_threshold_infec_optim_strat_vec_animallevel_indivpersp,...
                                    optim_strat_ternary_flag,...
                                    ternary_axes_labels,...
                                    title_string_indiv_persp,...
                                    strat_exceed_threshold_infections_chosen_cmap,...
                                    strat_exceed_threshold_infections_colourbar_label,...
                                    strat_exceed_threshold_infections_cbar_range,...
                                    strat_exceed_threshold_infections_colourbar_tick_pos,...
                                    strat_exceed_threshold_infections_colourbar_tick_labels,...
                                    label_fontsize,...
                                    tick_label_fontsize,...
                                    savefile_name_point_plot,...
                                    save_figure_flag)
        end
     end

    %%
    %--------------------------------------------------------------------------
    %%% Generate vaccinations under the optimal strategy plots for each relative cost of vacc
    %--------------------------------------------------------------------------
    if plot_infections_vaccs_flag == true
        optim_strat_ternary_flag =false;
        % Infection and vaccination stats not computed. Load from file instead
        if compute_infections_vacc_flag == false
            if strcmp(config_ID,'cumbria')
                load('OptimBehaviour_InfVacc_PlotData/cumbria_inf_vacc_ternary_plot_data.mat',...
                    'percentage_vacc_array_animallevel_popnpersp',...
                    'percentage_vacc_array_animallevel_indivpersp')
            elseif strcmp(config_ID,'cumbria_alt')
                load('OptimBehaviour_InfVacc_PlotData/cumbria_alt_inf_vacc_ternary_plot_data.mat',...
                    'percentage_vacc_array_animallevel_popnpersp',...
                    'percentage_vacc_array_animallevel_indivpersp')
            elseif strcmp(config_ID,'devon')
                load('OptimBehaviour_InfVacc_PlotData/devon_inf_vacc_ternary_plot_data.mat',...
                    'percentage_vacc_array_animallevel_popnpersp',...
                    'percentage_vacc_array_animallevel_indivpersp')
            elseif strcmp(config_ID,'devon_alt')
                load('OptimBehaviour_InfVacc_PlotData/devon_alt_inf_vacc_ternary_plot_data.mat',...
                    'percentage_vacc_array_animallevel_popnpersp',...
                    'percentage_vacc_array_animallevel_indivpersp')
            end
        end

        for rel_vacc_cost_itr = 1:length(rel_cost_of_vacc_idxs)

            % Assign index corresponding to relative vaccine cost in this iteration
            % to variable
            val_cost_of_vacc_col_idx = rel_cost_of_vacc_idxs(rel_vacc_cost_itr);

            % Set title for the plot
            rel_vacc_cost_val = (val_cost_of_vacc_col_idx*0.01)-0.01;
            % title_string = ['Relative cost of vaccination: ', num2str(rel_vacc_cost_val)];
            title_string_popn_persp = ['Population perspective, C_{V} = ', num2str(rel_vacc_cost_val)];
            title_string_indiv_persp = ['Individual perspective, C_{V} = ', num2str(rel_vacc_cost_val)];

            % Extract median vaccinations under the optimal strategy data
            median_vaccs_optim_strat_vec_animallevel_popnpersp = percentage_vacc_array_animallevel_popnpersp(val_cost_of_vacc_col_idx,:,2);
            median_vaccs_optim_strat_vec_animallevel_indivpersp = percentage_vacc_array_animallevel_indivpersp(val_cost_of_vacc_col_idx,:,2);

            %--  animallevel_popnpersp

            % Produce data point ternary plot
            savefile_name_point_plot = ['OptimBehaviour_VaccUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_popn_persp_sellke_run_median_vacc_ternary_data_point_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx),'.pdf'];
            data_point_ternary(ternay_plot_pos_array,...
                                    median_vaccs_optim_strat_vec_animallevel_popnpersp,...
                                    optim_strat_ternary_flag,...
                                    ternary_axes_labels,...
                                    title_string_popn_persp,...
                                    strat_vaccs_chosen_cmap,...
                                    strat_vaccs_colourbar_label,...
                                    strat_vaccs_cbar_range,...
                                    strat_vaccs_colourbar_tick_pos,...
                                    strat_vaccs_colourbar_tick_labels,...
                                    label_fontsize,...
                                    tick_label_fontsize,...
                                    savefile_name_point_plot,...
                                    save_figure_flag)

            %--  animallevel_indpersp

            % Produce data point ternary plot
            savefile_name_point_plot = ['OptimBehaviour_VaccUnderOptimStrat_SellkeRuns_TernaryPlots/animallevel_indiv_persp_sellke_run_median_vacc_ternary_data_point_plot_', config_ID ,'_vacc_cost_ID_', num2str(val_cost_of_vacc_col_idx),'.pdf'];
            data_point_ternary(ternay_plot_pos_array,...
                                    median_vaccs_optim_strat_vec_animallevel_indivpersp,...
                                    optim_strat_ternary_flag,...
                                    ternary_axes_labels,...
                                    title_string_indiv_persp,...
                                    strat_vaccs_chosen_cmap,...
                                    strat_vaccs_colourbar_label,...
                                    strat_vaccs_cbar_range,...
                                    strat_vaccs_colourbar_tick_pos,...
                                    strat_vaccs_colourbar_tick_labels,...
                                    label_fontsize,...
                                    tick_label_fontsize,...
                                    savefile_name_point_plot,...
                                    save_figure_flag)
        end
    end
end

%% SUPPORTING FUNCTIONS TO MAP CONTROL STRATEGY TO PLOT VALUE
function [data_array] = modify_value(data_array,original_val,plot_val)

    % For values of input_array matching no_control_original_val, replace with no_control_plot_val
    data_array(data_array==original_val) = plot_val;

end

%% SUPPORTING FUNCTION TO EXTRACT INFECTION & VACCINATION OUTCOMES
function [percentage_threshold_inf_exceeded_vec,...
            percentage_infected_vec,...
            percentage_vacc_vec] = get_infection_vacc_outcomes(rel_vacc_cost_idx,...
                                                                     vacc_strat_idx,...
                                                                     batch_ID_offset,...
                                                                     n_risk_measures,...
                                                                     epi_outputs_agg_directory,...
                                                                     optim_control_array,...
                                                                     prem_num,...
                                                                     threshold_infections_vals,...
                                                                     prctiles_infections_vaccs)
    % Premises level, population perspective
    optim_stratID = optim_control_array(rel_vacc_cost_idx,vacc_strat_idx) + 1;
        % optim_control_array(rel_vacc_cost_itr,vacc_strat_itr)
        % returns an integer value [0,1,...,n_risk_measures].
        % Add 1 to get the configuration within the scenario
        % set to be used.

    % Chosen batch ID
    chosen_batch_ID = batch_ID_offset + ((vacc_strat_idx-1)*n_risk_measures) + optim_stratID;

    % Load the data from the relevant Batch ID
    optim_strat_epi_outputs_filename = [epi_outputs_agg_directory,'PremPerDiseaseState_BatchID', num2str(chosen_batch_ID) ,'.txt'];
    optim_strat_epi_outputs_array = readmatrix(optim_strat_epi_outputs_filename);

    % For that batch ID, get the number of premises infected and vaccinated
    optim_strat_infections = optim_strat_epi_outputs_array(:,6);
    optim_strat_vaccs = optim_strat_epi_outputs_array(:,7);

    % Compute the number of replicates exceeding the specified infected
    % premises threshold values
    percentage_threshold_inf_exceeded_vec = (sum(optim_strat_infections >= threshold_infections_vals)/length(optim_strat_infections))*100;

    % Compute requested percentiles for infection data across all replicates
    % Turn into percentage based on total number of premises
    percentage_infected_vec = (prctile(optim_strat_infections,prctiles_infections_vaccs)/prem_num)*100;
    percentage_vacc_vec = (prctile(optim_strat_vaccs,prctiles_infections_vaccs)/prem_num)*100;
end

%% SUPPORTING FUNCTIONS TO CONSTRUCT TERNARY PLOT

%--------------------------------------------------------------------------
%%% Construct data point plot
%--------------------------------------------------------------------------
function data_point_ternary(ternay_plot_pos_array,...
                            optim_strat_vec,...
                            optim_strat_ternary_flag,...
                            ternary_axes_labels,...
                            title_string,...
                            chosen_cmap,...
                            colourbar_label,...
                            cbar_range,...
                            colourbar_tick_pos,...
                            colourbar_tick_labels,...
                            label_fontsize,...
                            tick_label_fontsize,...
                            save_filename,...
                            save_figure_flag)
    %Intialise figure
    position = [100, 100, 2*550, 1.1*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])

    % Plot the ternary axis system
    [h,hg,htick]=terplot;

    % Set the colormap
    colormap(chosen_cmap)

    % Plot the data. Call function based on whether plotting optimal
    % strategy data or a different summary statistic
    if optim_strat_ternary_flag == true
        [hd]=ternaryc_optim_strat(ternay_plot_pos_array(:,1),...
                            ternay_plot_pos_array(:,2),...
                            ternay_plot_pos_array(:,3),...
                            optim_strat_vec,...
                            'o');
    else
        [hd]=ternaryc_emh(ternay_plot_pos_array(:,1),...
                            ternay_plot_pos_array(:,2),...
                            ternay_plot_pos_array(:,3),...
                            optim_strat_vec,...
                            'o',...
                            cbar_range);
    end

    % Add the labels
    hlabels=terlabel(ternary_axes_labels{1},ternary_axes_labels{2},ternary_axes_labels{3});

    % Add the title
    title(title_string,...
          'Units', 'normalized',...
          'Position', [0.5, 1.0125, 0],...
          'FontSize',label_fontsize)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modify plot settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--  Change the color of the grid lines
    set(hg(:,3),'color','k')
    set(hg(:,2),'color','k')
    set(hg(:,1),'color','k')
    %--  Change the marker size
    set(hd,'markersize',15)
    %--  Modify the labels
    set(hlabels,'fontsize',label_fontsize)
    set(hlabels(3),'color','k')
    set(hlabels(2),'color','k')
    set(hlabels(1),'color','k')
    %--  Modify the tick labels
    htick_labels = ["0" "10" "20" "30" "40" "50" "60" "70" "80" "90" "100"];
    for label_itr = 1:11
        htick(label_itr,1).String = htick_labels(label_itr);
        htick(label_itr,2).String = htick_labels(label_itr);
        htick(label_itr,3).String = htick_labels(label_itr);
    end
    set(htick(:,1),'color','k','linewidth',3,'fontsize',tick_label_fontsize)
    set(htick(:,2),'color','k','linewidth',3,'fontsize',tick_label_fontsize)
    set(htick(:,3),'color','k','linewidth',3,'fontsize',tick_label_fontsize)

    %-- Move "100" tick label to left so does not spill onto marker
    set(htick(end,2),'Position',[-0.0879016858716176 -5.55111512312578e-17 0])

    %-- Move ticks on bottom axis to the left
    set(htick(1,1),'Position',get(htick(1,1),'Position') - [0.02 0 0])
    for tick_itr = 2:length(htick(:,1))
        set(htick(tick_itr,1),'Position',get(htick(tick_itr,1),'Position') - [0.03 0 0])
    end

    %--  Change the color of the patch
    set(h,'facecolor',[1. 1. 1.],'edgecolor','k')

    %-- Alter marker style (applicable to optimal strategy plot only)
    % Check if colour of marker handle is grey, corresponding to optimal
    % strategy of no vaccination.
    if optim_strat_ternary_flag == true
        for handle_itr = 1:length(hd)
            if isequal(hd(handle_itr).Color,[0.5 0.5 0.5])
                hd(handle_itr).Marker = 'x';
                hd(handle_itr).MarkerSize = 15;
                hd(handle_itr).LineWidth = 3;
            end
        end
    end

    %-- Add colourbar and label
    cc = colorbar;
    ylabel(cc, colourbar_label)
    % set(gca,'ColorScale','log')

    %-- Set colourbar limits
    caxis(cbar_range);

    %-- Set colourbar ticks
    cc.YTick = colourbar_tick_pos;
    cc.YTickLabel = colourbar_tick_labels;

    %-- Alter fontsize
    cc.FontSize = label_fontsize;

    %-- Set colourbar position
    cc_pos =  cc.Position; %gets the positon and size of the color bar
    set(cc,'Position',[cc_pos(1)+0.05 cc_pos(2) cc_pos(3) cc_pos(4)])% To change size

    %-- Modify plot properties
    set(gca,'LineWidth',1)
    box on

    % %--  Change the colorbar
    % set(hcb,'xcolor','w','ycolor','w')
    % %--  Modify the figure color
    % set(gcf,'color',[0 0 0.3])
    % %-- Change some defaults
    % set(gcf,'paperpositionmode','auto','inverthardcopy','off')

    % Save the figure
    if save_figure_flag == true
        % annotation('rectangle',[0.28 0.02 0.665 0.98],'Color','w'); % Add margin around the figure
        exportgraphics(gcf,save_filename,'BackgroundColor','none','ContentType','vector')
    end
end
