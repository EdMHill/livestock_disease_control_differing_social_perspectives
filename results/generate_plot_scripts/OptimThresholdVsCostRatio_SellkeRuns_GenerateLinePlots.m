% Purpose:
% Produce plot of optimal decision point against cost ratio
% 
% Outbreak simualtions used the Sellke construction.
%
% MATLAB version: R2021b
% Date: 3rd November 2021
%--------------------------------------------------------------------------

function OptimThresholdVsCostRatio_SellkeRuns_GenerateLinePlots(save_fig_flag,config_ID)

    %%
    %----------------------------------------------------------------------
    %%% Specify variables dependent on configuration in use 
    %%% Input filename specifier, save filename ID, plot titles & legend
    %%% inclusion
    %----------------------------------------------------------------------
    if strcmp(config_ID,'cumbria')
        input_filename_suffix = 'Cumbria_vacc_distance_risk_measure_scenID921';
        save_filename_suffix = 'OptimBehaviour_SellkeRuns_Plots/cumbria_scenID921_';
        county_specific_plot_title = 'Cumbria';
        include_legend_flag = false;
    elseif strcmp(config_ID,'cumbria_alt')
        input_filename_suffix = 'Cumbria_alt_pathogen_vacc_distance_risk_measure_scenID1221';
        save_filename_suffix = 'OptimBehaviour_SellkeRuns_Plots/cumbria_alt_scenID1221_';
        county_specific_plot_title = 'Cumbria';
        include_legend_flag = false;
    elseif strcmp(config_ID,'devon')
        input_filename_suffix = 'Devon_vacc_distance_risk_measure_scenID921';
        save_filename_suffix = 'OptimBehaviour_SellkeRuns_Plots/devon_scenID921_';
        county_specific_plot_title = 'Devon';
        include_legend_flag = true;
    elseif strcmp(config_ID,'devon_alt')
        input_filename_suffix = 'Devon_alt_pathogen_vacc_distance_risk_measure_scenID1221';
        save_filename_suffix = 'OptimBehaviour_SellkeRuns_Plots/devon_alt_scenID1221_';
        county_specific_plot_title = 'Devon';
        include_legend_flag = true;
    else
        error('Invalid config_ID provided.')
    end
    
    OptimThresholdVals_PremLevel_PopnPersp_FileName = ['OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_PremLevel_PopnPersp_' input_filename_suffix '.txt'];
    OptimThresholdVals_AnimalLevel_PopnPersp_FileName = ['OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_PopnPersp_' input_filename_suffix '.txt'];
    
    OptimThresholdVals_PremLevel_IndivPersp_FileName = ['OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_PremLevel_IndivPersp_' input_filename_suffix '.txt']';
    OptimThresholdVals_AnimalLevel_IndivPersp_FileName = ['OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_IndivPersp_' input_filename_suffix '.txt'];
    
    %----------------------------------------------------------------------
    %%% Load optimal strategy summary statistic data
    %----------------------------------------------------------------------
    OptimThresholdVals_PremLevel_PopnPersp = readmatrix(OptimThresholdVals_PremLevel_PopnPersp_FileName);
    OptimThresholdVals_AnimalLevel_PopnPersp = readmatrix(OptimThresholdVals_AnimalLevel_PopnPersp_FileName);
    
    OptimThresholdVals_PremLevel_IndivPersp = readmatrix(OptimThresholdVals_PremLevel_IndivPersp_FileName);
    OptimThresholdVals_AnimalLevel_IndivPersp = readmatrix(OptimThresholdVals_AnimalLevel_IndivPersp_FileName);
    
    %----------------------------------------------------------------------
    %%% Load infection and vaccination (under optimal strategy) summary statistic data
    %----------------------------------------------------------------------
    % Load the infection and vaccination data
    
    if strcmp(config_ID,'cumbria')
        infection_vacc_data = load('OptimBehaviour_InfVacc_PlotData/cumbria_inf_vacc_ternary_plot_data.mat');
    elseif strcmp(config_ID,'cumbria_alt')
        infection_vacc_data = load('OptimBehaviour_InfVacc_PlotData/cumbria_alt_inf_vacc_ternary_plot_data.mat');
    elseif strcmp(config_ID,'devon')
        infection_vacc_data = load('OptimBehaviour_InfVacc_PlotData/devon_inf_vacc_ternary_plot_data.mat');
    elseif strcmp(config_ID,'devon_alt')
        infection_vacc_data = load('OptimBehaviour_InfVacc_PlotData/devon_alt_inf_vacc_ternary_plot_data.mat');
    else
        error('Invalid config_ID provided.')
    end
    
    % Assign relevant variables from data structure to plotting variables
    % Column 21 for scenario 21 (all reactionary scenario)
        % Infection, slice to pick out the threshold count to match or exceed: threshold_infections_vals = [10 25 50 100 250 500 1000];
    infection_exceed_thereshold_popnpersp = infection_vacc_data.percentage_threshold_inf_exceeded_array_animallevel_popnpersp(:,21,2);
    infection_exceed_thereshold_indivpersp = infection_vacc_data.percentage_threshold_inf_exceeded_array_animallevel_indivpersp(:,21,2);
    
        % Percentage of premises infected, slice to pick out the percentile: [2.5 50 97.5] e.g. Index 2 for median.
    infection_data_popnpersp_median = infection_vacc_data.percentage_infected_array_animallevel_popnpersp(:,21,2);
    infection_data_indivpersp_median = infection_vacc_data.percentage_infected_array_animallevel_indivpersp(:,21,2);
    
    infection_data_popnpersp_LB = infection_vacc_data.percentage_infected_array_animallevel_popnpersp(:,21,1);
    infection_data_indivpersp_LB = infection_vacc_data.percentage_infected_array_animallevel_indivpersp(:,21,1);
    
    infection_data_popnpersp_UB = infection_vacc_data.percentage_infected_array_animallevel_popnpersp(:,21,3);
    infection_data_indivpersp_UB = infection_vacc_data.percentage_infected_array_animallevel_indivpersp(:,21,3);

        % Percentage of premises vaccinated, slice to pick out the percentile: [2.5 50 97.5] e.g. Index 2 for median.
    vacc_data_popnpersp_median = infection_vacc_data.percentage_vacc_array_animallevel_popnpersp(:,21,2);
    vacc_data_indivpersp_median = infection_vacc_data.percentage_vacc_array_animallevel_indivpersp(:,21,2);
    
    vacc_data_popnpersp_LB = infection_vacc_data.percentage_vacc_array_animallevel_popnpersp(:,21,1);
    vacc_data_indivpersp_LB = infection_vacc_data.percentage_vacc_array_animallevel_indivpersp(:,21,1);
    
    vacc_data_popnpersp_UB = infection_vacc_data.percentage_vacc_array_animallevel_popnpersp(:,21,3);
    vacc_data_indivpersp_UB = infection_vacc_data.percentage_vacc_array_animallevel_indivpersp(:,21,3);
    
    %----------------------------------------------------------------------
    %%% Declare vaccine to infection cost ratio values that were tested
    %----------------------------------------------------------------------
    VaccToInfCostRatio = 0:0.01:1; %We set C_I = 1
    
    %%
    %----------------------------------------------------------------------
    %%% Modify data if needed
    %----------------------------------------------------------------------
    
    % For values labelled as no control, map to value at bottom or top of
    % y-axis range as appropriate
    no_control_original_val = 0;
    no_control_plot_val = 0;
    
    OptimThresholdVals_PremLevel_PopnPersp = modify_value(OptimThresholdVals_PremLevel_PopnPersp,no_control_original_val,no_control_plot_val);
    OptimThresholdVals_AnimalLevel_PopnPersp = modify_value(OptimThresholdVals_AnimalLevel_PopnPersp,no_control_original_val,no_control_plot_val);
    
    OptimThresholdVals_PremLevel_IndivPersp = modify_value(OptimThresholdVals_PremLevel_IndivPersp,no_control_original_val,no_control_plot_val);
    OptimThresholdVals_AnimalLevel_IndivPersp = modify_value(OptimThresholdVals_AnimalLevel_IndivPersp,no_control_original_val,no_control_plot_val);
    
    % For values involving control, map original value to plot value
    control_original_vals = 1:1:10;
    control_plot_vals = 1:1:10;
    % original_vals = [0.5,1,3,5,10];
    % plot_vals = [1,2,3,4,5];
    for alter_val_idx = 1:length(control_original_vals)
        current_itr_original_val = control_original_vals(alter_val_idx); 
        current_itr_plot_val = control_plot_vals(alter_val_idx);
    
        OptimThresholdVals_PremLevel_PopnPersp = modify_value(OptimThresholdVals_PremLevel_PopnPersp,current_itr_original_val,current_itr_plot_val);
        OptimThresholdVals_AnimalLevel_PopnPersp = modify_value(OptimThresholdVals_AnimalLevel_PopnPersp,current_itr_original_val,current_itr_plot_val);
        
        OptimThresholdVals_PremLevel_IndivPersp = modify_value(OptimThresholdVals_PremLevel_IndivPersp,current_itr_original_val,current_itr_plot_val);
        OptimThresholdVals_AnimalLevel_IndivPersp = modify_value(OptimThresholdVals_AnimalLevel_IndivPersp,current_itr_original_val,current_itr_plot_val);
    end
    
    %%
    %----------------------------------------------------------------------
    %%% Set up plotting variables
    %----------------------------------------------------------------------
    
    % Axes labels
    yaxis_name = 'Risk threshold: Notified premises distance (km)';
    yaxis_name_infection = 'Replicates with 25+ premises infected (%)';
    yaxis_name_infection_with_PI = 'Premises infected (%)';
    yaxis_name_vacc = 'Premises vaccinated (%)';
    
    xaxis_name = 'Relative cost of vaccination, C_V';
    
    % y-axis limits & tick labels - optimal strategy plots
    yaxis_lims = [0 10];
    yaxis_tick_pos = 0:1:10;
    yaxis_tick_labels = {'0','1','2','3','4','5','6','7','8','9','10'};
    %yaxis_tick_pos = [0 1 2 3 4 5];
    %yaxis_tick_labels = {'0.0','0.5','1.0','3.0','5.0','10.0'};
    
    % y-axis limits & tick labels - infection plots (median profile only)
    yaxis_lims_infection = [0 100];
    yaxis_tick_pos_infection = 0:10:100;
    yaxis_tick_labels_infection = 0:10:100;

    % y-axis limits & tick labels - infection plots (with prediction intervals)
    yaxis_lims_infection_with_PI = [0 42];
    yaxis_tick_pos_infection_with_PI = 0:4:40;
    yaxis_tick_labels_infection_with_PI = 0:4:40;
    
    % y-axis limits & tick labels - vaccination plots (median profile only)
    yaxis_lims_vacc = [0 10];
    yaxis_tick_pos_vacc = 0:1:10;
    yaxis_tick_labels_vacc = 0:1:10;
    
    % y-axis limits & tick labels - vaccination plots (with prediction intervals)
    yaxis_lims_vacc_with_PI = [0 20];
    yaxis_tick_pos_vacc_with_PI = 0:2:20;
    yaxis_tick_labels_vacc_with_PI = 0:2:20;
    
    % x-axis tick labels
    xaxis_tick_pos = 0:0.1:1;
    xaxis_tick_labels = 0:0.1:1;
    
    % Set plot face transparency
    facealpha_val = 0.2;

    % Fontsize
    label_fontsize = 22;
    
    %%
    %----------------------------------------------------------------------
    %%% Assign select data to variables
    %----------------------------------------------------------------------
    
    %Note, summary statistics output are the mean and the following quantiles
    % [0.25,0.5,0.975,0,1,0.25,0.75,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9]
    
    
    %Pick out median values
    % OptimThresholdVals_PremLevel_PopnPersp_Median = OptimThresholdVals_PremLevel_PopnPersp(:,3);
    OptimThresholdVals_AnimalLevel_PopnPersp_Median = OptimThresholdVals_AnimalLevel_PopnPersp(:,3);
    
    % OptimThresholdVals_PremLevel_IndivPersp_Median = OptimThresholdVals_PremLevel_IndivPersp(:,3);
    OptimThresholdVals_AnimalLevel_IndivPersp_Median = OptimThresholdVals_AnimalLevel_IndivPersp(:,3);
    
%     
%     %Get lower bound and upper bounds
%     OptimThresholdVals_PremLevel_PopnPersp_95PI_LB = OptimThresholdVals_PremLevel_PopnPersp(:,2);
%     OptimThresholdVals_PremLevel_PopnPersp_95PI_UB = OptimThresholdVals_PremLevel_PopnPersp(:,4);
%     
%     OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB = OptimThresholdVals_AnimalLevel_PopnPersp(:,2);
%     OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB = OptimThresholdVals_AnimalLevel_PopnPersp(:,4);
%     
%     OptimThresholdVals_PremLevel_IndivPersp_95PI_LB = OptimThresholdVals_PremLevel_IndivPersp(:,2);
%     OptimThresholdVals_PremLevel_IndivPersp_95PI_UB = OptimThresholdVals_PremLevel_IndivPersp(:,4);
%     
%     OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB = OptimThresholdVals_AnimalLevel_IndivPersp(:,2);
%     OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB = OptimThresholdVals_AnimalLevel_IndivPersp(:,4);
%     
%     %Get minimum and maximum
%     OptimThresholdVals_PremLevel_PopnPersp_min = OptimThresholdVals_PremLevel_PopnPersp(:,5);
%     OptimThresholdVals_PremLevel_PopnPersp_max = OptimThresholdVals_PremLevel_PopnPersp(:,6);
%     
%     OptimThresholdVals_AnimalLevel_PopnPersp_min = OptimThresholdVals_AnimalLevel_PopnPersp(:,5);
%     OptimThresholdVals_AnimalLevel_PopnPersp_max = OptimThresholdVals_AnimalLevel_PopnPersp(:,6);
%     
%     OptimThresholdVals_PremLevel_IndivPersp_min = OptimThresholdVals_PremLevel_IndivPersp(:,5);
%     OptimThresholdVals_PremLevel_IndivPersp_max = OptimThresholdVals_PremLevel_IndivPersp(:,6);
%     
%     OptimThresholdVals_AnimalLevel_IndivPersp_min = OptimThresholdVals_AnimalLevel_IndivPersp(:,5);
%     OptimThresholdVals_AnimalLevel_IndivPersp_max = OptimThresholdVals_AnimalLevel_IndivPersp(:,6);
    
    %%
    %----------------------------------------------------------------------
    %%% Construct optimal control line profiles (displays costs applied at animal level only)
    %----------------------------------------------------------------------
    
    %Intialise figure
    position = [100, 100, 1.2*550, 1.2*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(); 
    clf;
    set(fig,'Color', [1 1 1])
    hold on
    
    %Add population perspective data.
    p1 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_Median,...
                    'Color',[0 0 0.8],...
                    'LineWidth',1.5,...
                    'DisplayName','Population perspective');
       
    %Add farmer perspective data.
    p2 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_Median,...
                    'Color',[0.8 0 0],...
                    'LineStyle','--',...
                    'LineWidth',1.5,...
                    'DisplayName','Individual perspective');
    
    %Title
    title(county_specific_plot_title)
    
    %X-axis properties
    xlabel(xaxis_name)
    xticks(xaxis_tick_pos)
    xticklabels(xaxis_tick_labels)
    
    %Y-axis properties
    ylabel(yaxis_name)
    ylim(yaxis_lims)
    yticks(yaxis_tick_pos)
    yticklabels(yaxis_tick_labels)
    
    %Add legend
    if include_legend_flag == true
        legend([p1;p2],'Location','Northeast',...
                        'Fontsize',label_fontsize)
    end
    
    %Axes properties
    set(gca,'FontSize',label_fontsize)
    set(gca,'LineWidth',1)
    box on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Save figure to file
    FileName = [save_filename_suffix,'optimise_median_cost_animallevel_sellke.pdf'];
    if save_fig_flag == true
        exportgraphics(gcf,FileName,'BackgroundColor','none','ContentType','vector')
        % export_fig(FileName,'-pdf','-transparent','-painters','-r1200')
    end
    
    %%
    %----------------------------------------------------------------------
    %%% Construct infection & vaccination line profiles under optimal control strategy 
    %%% (displays costs applied at animal level only)
    %----------------------------------------------------------------------
    
    %%% INFECTION FIGURE - START %%%
    
    %Intialise figure
    position = [100, 100, 1.2*550, 1.2*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(); 
    clf;
    set(fig,'Color', [1 1 1])
    hold on
    
    %Add population perspective data.
    plot(VaccToInfCostRatio,infection_exceed_thereshold_popnpersp,...
                    'Color',[0 0 0.8],...
                    'LineWidth',1.5,...
                    'DisplayName','Population perspective');
       
    %Add farmer perspective data.
    plot(VaccToInfCostRatio,infection_exceed_thereshold_indivpersp,...
                    'Color',[0.8 0 0],...
                    'LineStyle','--',...
                    'LineWidth',1.5,...
                    'DisplayName','Individual perspective');
    
    %Title
    title(county_specific_plot_title)
    
    %X-axis properties
    xlabel(xaxis_name)
    xticks(xaxis_tick_pos)
    xticklabels(xaxis_tick_labels)
    
    %Y-axis properties
    ylabel(yaxis_name_infection)
    ylim(yaxis_lims_infection)
    yticks(yaxis_tick_pos_infection)
    yticklabels(yaxis_tick_labels_infection)
    
    % %Add legend
    % if include_legend_flag == true
    %     legend([p1;p2],'Location','Northeast',...
    %                     'Fontsize',20)
    % end
    
    %Axes properties
    set(gca,'FontSize',label_fontsize)
    set(gca,'LineWidth',1)
    box on
    
    % Save figure to file
    FileName = [save_filename_suffix,'all_reactionary_infection_burden_profile_sellke.pdf'];
    if save_fig_flag == true
        exportgraphics(gcf,FileName,'BackgroundColor','none','ContentType','vector')
        % export_fig(FileName,'-pdf','-transparent','-painters','-r1200')
    end
    
    %%% INFECTION FIGURE - END %%%
    
    %%% VACC FIGURE - START %%%
    
    %Intialise figure
    position = [100, 100, 1.2*550, 1.2*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(); 
    clf;
    set(fig,'Color', [1 1 1])
    hold on
    
    %Add population perspective data.
    plot(VaccToInfCostRatio,vacc_data_popnpersp_median,...
                    'Color',[0 0 0.8],...
                    'LineWidth',1.5,...
                    'DisplayName','Population perspective');
       
    %Add farmer perspective data.
    plot(VaccToInfCostRatio,vacc_data_indivpersp_median,...
                    'Color',[0.8 0 0],...
                    'LineStyle','--',...
                    'LineWidth',1.5,...
                    'DisplayName','Individual perspective');
    
    %Title
    title(county_specific_plot_title)
    
    %X-axis properties
    xlabel(xaxis_name)
    xticks(xaxis_tick_pos)
    xticklabels(xaxis_tick_labels)
    
    %Y-axis properties
    ylabel(yaxis_name_vacc)
    ylim(yaxis_lims_vacc)
    yticks(yaxis_tick_pos_vacc)
    yticklabels(yaxis_tick_labels_vacc)
    
%     %Add legend
%     if include_legend_flag == true
%         legend([p1;p2],'Location','Northeast',...
%                         'Fontsize',20)
%     end
    
    %Axes properties
    set(gca,'FontSize',label_fontsize)
    set(gca,'LineWidth',1)
    box on
    
    % Save figure to file
    FileName = [save_filename_suffix,'all_reactionary_vacc_profile_sellke.pdf'];
    if save_fig_flag == true
        exportgraphics(gcf,FileName,'BackgroundColor','none','ContentType','vector')
        %export_fig(FileName,'-pdf','-transparent','-painters','-r1200')
    end
    %%% VACC FIGURE - END %%%
    
    %%
    %----------------------------------------------------------------------
    %%% Construct infection line profiles under optimal control strategy with
    %%% prediction intervals
    %%% (displays costs applied at animal level only)
    %----------------------------------------------------------------------
    
    % To shade region between two curves, use the format
    % plot(x,y1)
    % plot(x,y2)
    % fill([x fliplr(x)], [y1 fliplr(y2)], 'g')
    % y2 has to be a row vector for fliplr to have an effect.
    
    % Set up x values for use in fill command 
    x2 = [VaccToInfCostRatio fliplr(VaccToInfCostRatio)];
    
    % Define variables giving boundaries of region to be filled
    inBetween_Infections_AnimalLevel_PopnPersp = [infection_data_popnpersp_UB' fliplr(infection_data_popnpersp_LB')];
    inBetween_Infections_AnimalLevel_IndivPersp = [infection_data_indivpersp_UB' fliplr(infection_data_indivpersp_LB')];
    
    %Intialise figure
    position = [100, 100, 1.2*550, 1.2*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(); 
    clf;
    set(fig,'Color', [1 1 1])
    hold on
    
    %Add population perspective data.
    fill(x2, inBetween_Infections_AnimalLevel_PopnPersp, [0 0 0.8],...
                 'FaceAlpha',facealpha_val);
    plot(VaccToInfCostRatio,infection_data_popnpersp_median,...
                    'Color',[0 0 0.8],...
                    'LineWidth',1.5,...
                    'DisplayName','Population perspective');
       
    %Add farmer perspective data.
    fill(x2, inBetween_Infections_AnimalLevel_IndivPersp, [0.8 0 0],...
                 'FaceAlpha',0.2);
    plot(VaccToInfCostRatio,infection_data_indivpersp_median,...
                    'Color',[0.8 0 0],...
                    'LineStyle','--',...
                    'LineWidth',1.5,...
                    'DisplayName','Individual perspective');
    
    %Title
    title(county_specific_plot_title)
    
    %X-axis properties
    xlabel(xaxis_name)
    xticks(xaxis_tick_pos)
    xticklabels(xaxis_tick_labels)
    
    %Y-axis properties
    ylabel(yaxis_name_infection_with_PI)
    ylim(yaxis_lims_infection_with_PI)
    yticks(yaxis_tick_pos_infection_with_PI)
    yticklabels(yaxis_tick_labels_infection_with_PI)
    
    %Axes properties
    set(gca,'FontSize',label_fontsize)
    set(gca,'LineWidth',1)
    box on
    
    % Save figure to file
    FileName = [save_filename_suffix,'all_reactionary_infection_profile_with_PI_sellke.pdf'];
    if save_fig_flag == true
        annotation('rectangle',[0.04 0 0.88 0.98],'Color','w'); % Add margin around the figure
        exportgraphics(gcf,FileName,'BackgroundColor','none','ContentType','vector')
    end

    %%
    %----------------------------------------------------------------------
    %%% Construct vaccination line profiles under optimal control strategy with
    %%% prediction intervalss
    %%% (displays costs applied at animal level only)
    %----------------------------------------------------------------------
    
    % To shade region between two curves, use the format
    % plot(x,y1)
    % plot(x,y2)
    % fill([x fliplr(x)], [y1 fliplr(y2)], 'g')
    % y2 has to be a row vector for fliplr to have an effect.
    
    % Set up x values for use in fill command 
    x2 = [VaccToInfCostRatio fliplr(VaccToInfCostRatio)];
    
    % Define variables giving boundaries of region to be filled
    inBetween_Vacc_AnimalLevel_PopnPersp = [vacc_data_popnpersp_UB' fliplr(vacc_data_popnpersp_LB')];
    inBetween_Vacc_AnimalLevel_IndivPersp = [vacc_data_indivpersp_UB' fliplr(vacc_data_indivpersp_LB')];
    
    %Intialise figure
    position = [100, 100, 1.2*550, 1.2*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(); 
    clf;
    set(fig,'Color', [1 1 1])
    hold on
    
    %Add population perspective data.
    fill(x2, inBetween_Vacc_AnimalLevel_PopnPersp, [0 0 0.8],...
                 'FaceAlpha',facealpha_val);
    plot(VaccToInfCostRatio,vacc_data_popnpersp_median,...
                    'Color',[0 0 0.8],...
                    'LineWidth',1.5,...
                    'DisplayName','Population perspective');
       
    %Add farmer perspective data.
    fill(x2, inBetween_Vacc_AnimalLevel_IndivPersp, [0.8 0 0],...
                 'FaceAlpha',0.2);
    plot(VaccToInfCostRatio,vacc_data_indivpersp_median,...
                    'Color',[0.8 0 0],...
                    'LineStyle','--',...
                    'LineWidth',1.5,...
                    'DisplayName','Individual perspective');
    
    %Title
    title(county_specific_plot_title)
    
    %X-axis properties
    xlabel(xaxis_name)
    xticks(xaxis_tick_pos)
    xticklabels(xaxis_tick_labels)
    
    %Y-axis properties
    ylabel(yaxis_name_vacc)
    ylim(yaxis_lims_vacc_with_PI)
    yticks(yaxis_tick_pos_vacc_with_PI)
    yticklabels(yaxis_tick_labels_vacc_with_PI)
    
    %Axes properties
    set(gca,'FontSize',label_fontsize)
    set(gca,'LineWidth',1)
    box on
    
    % Save figure to file
    FileName = [save_filename_suffix,'all_reactionary_vacc_profile_with_PI_sellke.pdf'];
    if save_fig_flag == true
        annotation('rectangle',[0.04 0 0.88 0.98],'Color','w'); % Add margin around the figure
        exportgraphics(gcf,FileName,'BackgroundColor','none','ContentType','vector')
    end
end

%%
% %%
% %----------------------------------------------------------------------
% %%% Construct plots (displays costs applied at premises and animal level)
% %----------------------------------------------------------------------
% 
% %To shade region between two curves, use the format
% %plot(x,y1)
% %plot(x,y2)
% %fill([x fliplr(x)], [y1 fliplr(y2)], 'g')
% % y2 has to be a row vector for fliplr to have an effect.
% 
% % Set up x values for use in fill command 
% x2 = [VaccToInfCostRatio fliplr(VaccToInfCostRatio)];
% 
% % Set up regions to fill
% inBetween_PremLevel_PopnPersp = [OptimThresholdVals_PremLevel_PopnPersp_max' fliplr(OptimThresholdVals_PremLevel_PopnPersp_min')];
% inBetween_PremLevel_IndPersp = [OptimThresholdVals_PremLevel_IndivPersp_max' fliplr(OptimThresholdVals_PremLevel_IndivPersp_min')];
% 
% inBetween_AnimalLevel_PopnPersp = [OptimThresholdVals_AnimalLevel_PopnPersp_max' fliplr(OptimThresholdVals_AnimalLevel_PopnPersp_min')];
% inBetween_AnimalLevel_IndPersp = [OptimThresholdVals_AnimalLevel_IndivPersp_max' fliplr(OptimThresholdVals_AnimalLevel_IndivPersp_min')];

% inBetween_PremLevel_PopnPersp = [OptimThresholdVals_PremLevel_PopnPersp_95PI_UB' fliplr(OptimThresholdVals_PremLevel_PopnPersp_95PI_LB')];
% inBetween_PremLevel_IndPersp = [OptimThresholdVals_PremLevel_IndivPersp_95PI_UB' fliplr(OptimThresholdVals_PremLevel_IndivPersp_95PI_LB')];
% 
% inBetween_AnimalLevel_PopnPersp = [OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB' fliplr(OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB')];
% inBetween_AnimalLevel_IndPersp = [OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB' fliplr(OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB')];

% %%
% %%% Dual panel (both perspectives), with 95% interval %%%
% 
% %Intialise figure
% position = [100, 100, 2*550, 450];
% set(0, 'DefaultFigurePosition', position);
% fig = figure(); 
% clf;
% set(fig,'Color', [1 1 1])
% 
% %%% Premises-level panel %%%
% subplot(1,2,1)
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_PremLevel_PopnPersp_95PI_UB;fliplr(OptimThresholdVals_PremLevel_PopnPersp_95PI_LB)],...
% %              'c',...
% %             'FaceAlpha',0.2)
% fill(x2, inBetween_PremLevel_PopnPersp, [0 0 0.8],...
%             'FaceAlpha',facealpha_val);
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_Median,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective')
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)  
% 
% 
% % 
% % % %Add farmer perspective data. Shade prediction interval
% fill(x2, inBetween_PremLevel_IndPersp, [0.8 0 0],...
%             'FaceAlpha',0.2);
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_Median,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective')
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)  
% 
% 
% %Title
% title('Premises-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% %Add legend
% legend([p1;p2],'Location','Northeast',...
%                 'Fontsize',18)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Animal-level panel %%% 
% subplot(1,2,2)
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% fill(x2, inBetween_AnimalLevel_PopnPersp, [0 0 0.8],...
%             'FaceAlpha',0.2);  
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_Median,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective');
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)  
% 
% 
%    
% %Add farmer perspective data. Shade prediction interval
% fill(x2, inBetween_AnimalLevel_IndPersp, [0.8 0 0],...
%             'FaceAlpha',0.2); 
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_Median,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective');
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)  
% 
% %Title
% title('Animal-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Save figure to file
% FileName = 'OptimBehaviour_Plots/CumbriaExamplePlot_optimise_cost_with_prediction_interval';
% if save_fig_flag == true
%     %export_fig(FileName,'-pdf','-transparent','-painters','-r1200')
% end
%%
% %%% Dual panel (both perspectives), optimise based on median cost %%%
% 
% %Intialise figure
% position = [100, 100, 2*550, 450];
% set(0, 'DefaultFigurePosition', position);
% fig = figure(); 
% clf;
% set(fig,'Color', [1 1 1])
% 
% %%% Premises-level panel %%%
% subplot(1,2,1)
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_PremLevel_PopnPersp_95PI_LB;fliplr(OptimThresholdVals_PremLevel_PopnPersp_95PI_UB)],...
% %              'c',...
% %             'FaceAlpha',0.2)
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_Median,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective')
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)  
% 
% 
% 
% % %Add farmer perspective data. Shade prediction interval
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_PremLevel_IndivPersp_95PI_LB;fliplr(OptimThresholdVals_PremLevel_IndivPersp_95PI_UB)],...
% %             'g')
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_Median,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective')
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)  
% % 
% 
% 
% %Title
% title('Premises-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% %Add legend
% legend([p1;p2],'Location','Northeast',...
%                 'Fontsize',18)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Animal-level panel %%% 
% subplot(1,2,2)
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB;fliplr(OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB)],...
% %             'c',...
% %             'FaceAlpha',0.2)   
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_Median,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective');
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)  
% 
% 
%    
% %Add farmer perspective data. Shade prediction interval
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB;fliplr(OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB)],...
% %             'g')
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_Median,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective');
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)  
% 
% %Title
% title('Animal-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Save figure to file
% FileName = [save_filename_suffix,'optimise_median_cost_sellke'];
% if save_fig_flag == true
%     export_fig(FileName,'-pdf','-transparent','-painters','-r1200')
% end
% %%
% %%% Dual panel (both perspectives), optimise based on maximum %%%
% 
% %Intialise figure
% position = [100, 100, 2*550, 450];
% set(0, 'DefaultFigurePosition', position);
% fig = figure(); 
% clf;
% set(fig,'Color', [1 1 1])
% 
% %%% Premises-level panel %%%
% subplot(1,2,1)
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_max,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective')
% 
% % %Add farmer perspective data. Shade prediction interval
% plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_max,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective')
% 
% %Title
% title('Premises-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Animal-level panel %%% 
% subplot(1,2,2)
% hold on
% 
% %Add national perspective data.
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_max,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective');
%    
% %Add farmer perspective data
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_max,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective');
% 
% 
% %Title
% title('Animal-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% % %Add legend
% % legend([p1;p2],'Location','Northwest',...
% %                 'Fontsize',18)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Save figure to file
% FileName = [save_filename_suffix,'optimise_maximum_sellke'];
% if save_fig_flag == true
%     export_fig(FileName,'-pdf','-transparent','-painters','-r1200')
% end
% %%
% %%% Dual panel (both perspectives), optimise based on minimum %%%
% 
% %Intialise figure
% position = [100, 100, 2*550, 450];
% set(0, 'DefaultFigurePosition', position);
% fig = figure(); 
% clf;
% set(fig,'Color', [1 1 1])
% 
% %%% Premises-level panel %%%
% subplot(1,2,1)
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_min,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective')
% 
% % %Add farmer perspective data. Shade prediction interval
% plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_min,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective')
% 
% %Title
% title('Premises-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Animal-level panel %%% 
% subplot(1,2,2)
% hold on
% 
% %Add national perspective data.
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_min,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective');
%    
% %Add farmer perspective data
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_min,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective');
% 
% 
% %Title
% title('Animal-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% % %Add legend
% % legend([p1;p2],'Location','Northwest',...
% %                 'Fontsize',18)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Save figure to file
% FileName = [save_filename_suffix,'optimise_minimum_sellke'];
% if save_fig_flag == true
%     export_fig(FileName,'-pdf','-transparent','-painters','-r1200')
% end

% %%
% %%% Dual panel (both perspectives), optimise based on 97.5th percentile %%%
% 
% %Intialise figure
% position = [100, 100, 2*550, 450];
% set(0, 'DefaultFigurePosition', position);
% fig = figure(); 
% clf;
% set(fig,'Color', [1 1 1])
% 
% %%% Premises-level panel %%%
% subplot(1,2,1)
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_95PI_UB,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective')
% 
% % %Add farmer perspective data. Shade prediction interval
% plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_95PI_UB,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective')
% 
% %Title
% title('Premises-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Animal-level panel %%% 
% subplot(1,2,2)
% hold on
% 
% %Add national perspective data.
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective');
%    
% %Add farmer perspective data
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective');
% 
% 
% %Title
% title('Animal-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% % %Add legend
% % legend([p1;p2],'Location','Northwest',...
% %                 'Fontsize',18)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Save figure to file
% FileName = 'OptimBehaviour_Plots/CumbriaExamplePlot_optimise_UB_percentile';
% export_fig(FileName,'-pdf','-transparent','-painters','-r1200')
% 
% %%
% %%% Dual panel (both perspectives), optimise based on 2.5th percentile %%%
% 
% %Intialise figure
% position = [100, 100, 2*550, 450];
% set(0, 'DefaultFigurePosition', position);
% fig = figure(); 
% clf;
% set(fig,'Color', [1 1 1])
% 
% %%% Premises-level panel %%%
% subplot(1,2,1)
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_95PI_LB,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective')
% 
% % %Add farmer perspective data. Shade prediction interval
% plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_95PI_LB,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective')
% 
% %Title
% title('Premises-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Animal-level panel %%% 
% subplot(1,2,2)
% hold on
% 
% %Add national perspective data.
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective');
%    
% %Add farmer perspective data
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective');
% 
% 
% %Title
% title('Animal-level costs')
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% % %Add legend
% % legend([p1;p2],'Location','Northwest',...
% %                 'Fontsize',18)
% 
% %Axes properties
% set(gca,'FontSize',18)
% set(gca,'LineWidth',1)
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Save figure to file
% FileName = 'OptimBehaviour_Plots/CumbriaExamplePlot_optimise_LB_percentile';
% export_fig(FileName,'-pdf','-transparent','-painters','-r1200')

% %%
% %%% Premises-level %%%
% 
% %Intialise figure
% position = [100, 100, 550, 450];
% set(0, 'DefaultFigurePosition', position);
% fig = figure(); 
% clf;
% set(fig,'Color', [1 1 1])
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_PremLevel_PopnPersp_95PI_LB;fliplr(OptimThresholdVals_PremLevel_PopnPersp_95PI_UB)],... 
% %             'c',...
% %             'FaceAlpha',0.2)
%         
% p1 = plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_Median,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5,...
%                 'DisplayName','Population perspective');
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_PopnPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)  
% 
%                 
%    
% %Add farmer perspective data. Shade prediction interval
% 
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_PremLevel_IndivPersp_95PI_LB;fliplr(OptimThresholdVals_PremLevel_IndivPersp_95PI_UB)],...
% %             'g',...
% %             'FaceAlpha',0.2)
% 
% p2 = plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_Median,...
%                 'Color',[0.8 0 0],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5,...
%                 'DisplayName','Individual perspective');
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',2.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_PremLevel_IndivPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',2.5)  
% 
% 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% %Add legend
% legend([p1;p2],'Location','Southwest')
% 
% %Axes properties
% set(gca,'FontSize',16)
% set(gca,'LineWidth',1)
% box on
% %%
% %%% Animal-level %%%
% 
% %Intialise figure
% position = [100, 100, 550, 450];
% set(0, 'DefaultFigurePosition', position);
% fig = figure(); 
% clf;
% set(fig,'Color', [1 1 1])
% hold on
% 
% %Add national perspective data. Shade prediction interval.
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB;fliplr(OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB)],...
% %             'c',...
% %             'FaceAlpha',0.2)
%                 
% plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_Median,...
%                 'Color',[0 0 0.8],...
%                 'LineWidth',1.5)
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_PopnPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineWidth',0.5)  
% 
% 
%    
%    
% % %Add farmer perspective data. Shade prediction interval
% % patch([VaccToInfCostRatio fliplr(VaccToInfCostRatio)],...
% %             [OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB;fliplr(OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB)],...
% %             'g')
% plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_Median,...
%                 'Color',[0 0 0.8],...
%                 'LineStyle','--',...
%                 'LineWidth',1.5)
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_95PI_LB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)        
% % plot(VaccToInfCostRatio,OptimThresholdVals_AnimalLevel_IndivPersp_95PI_UB,...
% %                     'Color',[0.5 0.5 0.5],...
% %                     'LineStyle','--',...
% %                     'LineWidth',0.5)  
% % 
% 
%                 
% %X-axis properties
% xlabel(xaxis_name)
% 
% %Y-axis properties
% ylabel(yaxis_name)
% ylim(yaxis_lims)
% yticks(yaxis_tick_pos)
% yticklabels(yaxis_tick_labels)
% 
% %Axes properties
% set(gca,'FontSize',16)
% set(gca,'LineWidth',1)
% box on


%%
%----------------------------------------------------------------------
%%% Supporting functions
%----------------------------------------------------------------------
function [data_array] = modify_value(data_array,original_val,plot_val)

    % For values of input_array matching no_control_original_val, replace with no_control_plot_val
    data_array(data_array==original_val) = plot_val;

end