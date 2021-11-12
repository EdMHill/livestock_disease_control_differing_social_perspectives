% Purpose:
% Script to call functions producing figures
%
% MATLAB version: R2021b
% Date: 3rd November 2021
%--------------------------------------------------------------------------

clear variables

%% Optimal strategy ternary plots (Figures 3-4 & S3-S12)

% Add colourbrewer and ternary plot packages to path
addpath('../../src/matlab_packages/cbrewer/')
addpath('../../src/matlab_packages/ternary2/')

% Specify if generated figures should be saved to file
save_fig_flag = false;

% Specify the configuration (location & pathogen attribute set) in use
config_ID = 'cumbria';
% config_ID = 'cumbria_alt';
% config_ID = 'devon';
% config_ID = 'devon_alt';

% Set scenario IDs & batch_ID_offset to be used
if strcmp(config_ID,'cumbria') || strcmp(config_ID,'devon')
    scen_IDs = 901:1:1131;
    batch_ID_offset = 15000;
elseif strcmp(config_ID,'cumbria_alt') || strcmp(config_ID,'devon_alt')
    scen_IDs = 1201:1:1431;
    batch_ID_offset = 20000;
else
    error('Invalid config_ID provided.')
end

% Relative cost of vaccination to generate plots for
rel_cost_of_vacc_idxs = [2,21,41,61,81,101]; % e.g. 1 for 0, 11 for 0.1, 51 for 0.5, 101 for 1.

% Flag to specify if strategy outputs should be computed and plotted
compute_optim_strategy_flag = true;
plot_optim_strategy_flag = false;

% Flag to specify if cost outputs should be produced
plot_costs_flag = false;

% Specify if infection & vaccination outputs should be plotted and what
% percentile values + exceedence of threshold infection should be computed
plot_infections_vaccs_flag = false;
compute_infections_vacc_flag = false; % If false, values are loaded from an existing MAT file, selected based on config_ID
prctiles_infections_vaccs = [50 2.5 97.5 5 95 25 75];
threshold_infections_vals = [10 25 50 100 250 500 1000];

% Call the function to produce the ternary plot set
[percentage_threshold_inf_exceeded_array_animallevel_popnpersp,...
          percentage_threshold_inf_exceeded_array_animallevel_indivpersp,...
          percentage_infected_array_animallevel_popnpersp,...
          percentage_infected_array_animallevel_indivpersp,...
          percentage_vacc_array_animallevel_popnpersp,...
          percentage_vacc_array_animallevel_indivpersp] = OptimThresholdVsCostRatio_SellkeRuns_GenerateTernaryPlot(save_fig_flag,...
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
                                                            threshold_infections_vals);

% Check if should save outputs to MAT file
if compute_infections_vacc_flag == true

    % Set save filename based on configuration
    if strcmp(config_ID,'cumbria')
        save_file_MAT = 'OptimBehaviour_InfVacc_PlotData/cumbria_inf_vacc_ternary_plot_data.mat';
    elseif strcmp(config_ID,'cumbria_alt')
        save_file_MAT = 'OptimBehaviour_InfVacc_PlotData/cumbria_alt_inf_vacc_ternary_plot_data.mat';
    elseif strcmp(config_ID,'devon')
        save_file_MAT = 'OptimBehaviour_InfVacc_PlotData/devon_inf_vacc_ternary_plot_data.mat';
    elseif strcmp(config_ID,'devon_alt')
        save_file_MAT = 'OptimBehaviour_InfVacc_PlotData/devon_alt_inf_vacc_ternary_plot_data.mat';
    else
        error('Invalid config_ID provided.')
    end

    save(save_file_MAT,'percentage_threshold_inf_exceeded_array_animallevel_popnpersp',...
                          'percentage_threshold_inf_exceeded_array_animallevel_indivpersp',...
                          'percentage_infected_array_animallevel_popnpersp',...
                          'percentage_infected_array_animallevel_indivpersp',...
                          'percentage_vacc_array_animallevel_popnpersp',...
                          'percentage_vacc_array_animallevel_indivpersp')
end

%% Example ternary plot (Figure S2)

% Add ternary plot package to path
addpath('../../src/matlab_packages/ternary2/')

% Specify if generated figures should be saved to file
save_fig_flag = true;

% Produce a single ternary plot. Markers on the corners. 
% One example somewhere in middle with grid lines highlighted
make_explainer_ternary_plot(save_fig_flag)