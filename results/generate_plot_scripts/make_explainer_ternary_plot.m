%Purpose:
% An explainer plot on reading values from a ternary plot
%
% Cite Ternary Plots package as
% Ulrich Theune (2021). Ternary Plots (https://www.mathworks.com/matlabcentral/fileexchange/7210-ternary-plots), 
% MATLAB Central File Exchange. Retrieved October 4, 2021.
%
% MATLAB version: R2021b
% Date: 12th November 2021
%--------------------------------------------------------------------------

function make_explainer_ternary_plot(save_fig_flag)

    %Intialise figure
    position = [100, 100, 1.5*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])

    % Font sizes
    label_fontsize = 22;
    tick_label_fontsize = 20;

    % Set marker colours
    all_precautionary_colour_vec = [0. 0 0.8];
    all_reactionary_colour_vec = [0 0. 0.8];
    all_nonvacc_colour_vec = [0 0 0.8];
    equalish_group_mix_colour_vec = [0.8 0. 0.];
    nonequal_group_mix_colour_vec = [0.8 0. 0.];

    % Plot the ternary axis system
    [h,hg,htick]=terplot;
    
    % Produce a single ternary plot. Markers on the corners. 
    % One example somewhere in middle with grid lines highlighted

    % Add markers
    [hd_1]=ternaryc(1,0,0,all_precautionary_colour_vec,'o');
    [hd_2]=ternaryc(0,1,0,all_reactionary_colour_vec,'d');
    [hd_3]=ternaryc(0,0,1,all_nonvacc_colour_vec,'v');
    [hd_4]=ternaryc(0.3,0.3,0.4,equalish_group_mix_colour_vec,'p');
    [hd_5]=ternaryc(0.1,0.7,0.2,nonequal_group_mix_colour_vec,'h');

    % Add the labels
    ternary_axes_labels = {'Precautionary (%)','Reactionary (%)','Non-vaccinators (%)'};
    hlabels=terlabel(ternary_axes_labels{1},ternary_axes_labels{2},ternary_axes_labels{3});
    
%     % Add the title
%     title(title_string,...
%         'Units', 'normalized',...
%         'Position', [0.5, 1.0125, 0],...
%         'FontSize',label_fontsize)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modify marker settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--  Change the marker size and colour
    set(hd_1,'MarkerSize',20,'MarkerFaceColor',all_precautionary_colour_vec,'MarkerEdgeColor',all_precautionary_colour_vec)
    set(hd_2,'MarkerSize',20,'MarkerFaceColor',all_reactionary_colour_vec,'MarkerEdgeColor',all_reactionary_colour_vec)
    set(hd_3,'MarkerSize',20,'MarkerFaceColor',all_nonvacc_colour_vec,'MarkerEdgeColor',all_nonvacc_colour_vec)
    set(hd_4,'MarkerSize',30,'MarkerFaceColor',equalish_group_mix_colour_vec,'MarkerEdgeColor',equalish_group_mix_colour_vec)
    set(hd_5,'MarkerSize',30,'MarkerFaceColor',nonequal_group_mix_colour_vec,'MarkerEdgeColor',nonequal_group_mix_colour_vec)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modify grid line settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--  Modify the original grid line colour
    set(hg(:,3),'Color','k','LineWidth',0.5)
    set(hg(:,2),'Color','k','LineWidth',0.5)
    set(hg(:,1),'Color','k','LineWidth',0.5)

    %--  Add lines for the equal(ish) group sizes
    set(hg(3,1),'Color','k','LineStyle','-','LineWidth',3,'XData',[0.3 0.5],'YData',[0 0.3464])
    set(hg(3,2),'Color','k','LineStyle',':','LineWidth',3,'XData',[0.35 0.5],'YData',[0.6062 0.3464])
    set(hg(4,3),'Color','k','LineStyle','--','LineWidth',3,'XData',[0.5 0.8],'YData',[0.3464 0.3464])

    %--  Alter grid lines for the nonequal group sizes
    set(hg(7,1),'Color','k','LineStyle','-','LineWidth',3,'XData',[0.7 0.8],'YData',[0 0.1732])
    set(hg(1,2),'Color','k','LineStyle',':','LineWidth',3,'XData',[0.45 0.8],'YData',[0.7794 0.1732])
    set(hg(2,3),'Color','k','LineStyle','--','LineWidth',3,'XData',[0.8 0.9],'YData',[0.1732 0.1732])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modify plot settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %-- Modify plot properties
    set(gca,'LineWidth',1)
    box on
    
    % If applicable, save file
    if save_fig_flag == true
        example_ternary_plot_save_filename = 'example_ternary_plots/example_ternary_plot.pdf';
        exportgraphics(gcf,example_ternary_plot_save_filename,'BackgroundColor','none','ContentType','vector')
    end

end