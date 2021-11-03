% Purpose:
% Script to call functions producing figures
%
% MATLAB version: R2021b
% Date: 3rd November 2021
%--------------------------------------------------------------------------

clear variables

%%
%--------------------------------------------------------------------------
% Set current directory and add paths
%--------------------------------------------------------------------------
cd('/Users/edhill/Documents/GitHub/LivestockInfectionFarmerLedControl/Results/sellke_construction/OptimBehavSimns_GenericModelSellke')
addpath('../../../src/Matlab_packages/cbrewer/')

%% Cumbria & Devon livestock population demographics (Figure 1)

% Set working directory
cd('/Users/edhill/Documents/GitHub/LivestockInfectionFarmerLedControl/Data/DataVis/2020_livestock_popn_data_plots/county_comparison/');

% Specify if generated figures should be saved to file
save_fig_flag = false;

% Set county_IDs as a vector
county_IDs = [8,10];
% county_IDs = [0,0];
            % Cumbria: 8
            % Devon: 10
            % Dummy data: 0

% Call the function to produce the plot set
between_county_livestock_comparison_script(save_fig_flag,county_IDs)

%% All reactionary vaccinators scenario plots (Figure 2)

% Set working directory
cd('/Users/edhill/Documents/GitHub/LivestockInfectionFarmerLedControl/Results/sellke_construction/OptimBehavSimns_GenericModelSellke')

% Specify if generated figures should be saved to file
save_fig_flag = false;

% Specify the configuration (location & pathogen attribute set) in use
config_ID = 'cumbria';
% config_ID = 'cumbria_alt';
% config_ID = 'devon';
% config_ID = 'devon_alt';

% Call the function to produce the plot set
OptimThresholdVsCostRatio_SellkeRuns_GenerateLinePlots(save_fig_flag,config_ID)

%% Optimal strategy ternary plots (Figure 3)

% Set working directory
cd('/Users/edhill/Documents/GitHub/LivestockInfectionFarmerLedControl/Results/sellke_construction/OptimBehavSimns_GenericModelSellke')

% Add ternary plot package to path
addpath('../../../src/Matlab_packages/ternary2/')

% Specify if generated figures should be saved to file
save_fig_flag = true;

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
compute_optim_strategy_flag = false;
plot_optim_strategy_flag = false;

% Flag to specify if cost outputs should be produced
plot_costs_flag = false;

% Specify if infection & vaccination outputs should be plotted and what
% percentile values + exceedence of threshold infection should be computed
plot_infections_vaccs_flag = true;
compute_infections_vacc_flag = false; % If false, values are loaded from an existing MAT file, selected based on config_ID
prctiles_infections_vaccs = [2.5 50 97.5];
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
