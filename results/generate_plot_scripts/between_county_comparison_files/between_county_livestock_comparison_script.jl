#=
Purpose:
Between county comparison of livestock (cattle and sheep) composition

Uses 2020 datasets

Julia version: 1.6.3
Date: 3rd November 2021
=#


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

"""
between_county_livestock_comparison_plots(save_figs_flag::Bool,
                                            county_IDs::Array{Int64,1})

Between county comparison of livestock (cattle and sheep) composition.

Inputs:
- `save_fig_flag::Bool`: Flag determining whether generated figures are saved to file or not.
- `county_IDs::Array{Int64,1}`: The county information to be analysed.

Outputs: None \n
Location: between\\_county\\_livestock\\_comparison\\_script.jl
"""
function between_county_livestock_comparison_plots(save_figs_flag::Bool,
                                                    county_IDs::Array{Int64,1})

    #----------------------------------------------------------------------
    # Set path to directory this file resides in
    #----------------------------------------------------------------------
    cd(dirname(@__FILE__))

    #----------------------------------------------------------------------
    # Specify number of counties to be analysed
    #----------------------------------------------------------------------
    n_counties = length(county_IDs)
        # county_IDs can be a vector
        # e.g. county_IDs = [8,10]; % Read in from vector
            # Cumbria: 8
            # Devon: 10
            # Dummy data: 0

    #----------------------------------------------------------------------
    # Specify global plot properties
    #----------------------------------------------------------------------

    # Group labels
    absolute_counts_groups = ["Cattle only","Sheep only","Cattle and sheep"]
    propns_groups = ["Cumbria","Devon"]

    # Legend labels
    absolute_counts_leg_labels = ["Cumbria" "Devon"]
    ecdf_leg_labels = ["Cumbria: cattle herds","Cumbria: sheep flocks","Devon: cattle herds","Devon: sheep flocks"]

    # Fontsize
    plot_fontsize_default = 14
    plot_fontsize_ecdfs = 12
    plot_fontsize_ecdfs_leg = 10
    plot_fontsize_map = 12

    # Set figure width and height
    fig_width = 1.2*550
    fig_height = 1.2*450

    fig_width_ecdfs = 1.0*550
    fig_height_ecdfs = 1.0*450

    # Colourmaps
    absolute_counts_colourmap = [RGB(0.2,0.2,0.2) RGB(0.8,0.8,0.8)]

    propns_colourmap = [ColorSchemes.Oranges_3.colors[1] ColorSchemes.Oranges_3.colors[2] ColorSchemes.Oranges_3.colors[3]]

    ecdfs_colourmap = [RGB(0.,0.,0.8) RGB(0.8,0.,0.)]

    # marker style
    ecdfs_marker_type = [:circle,:dtriangle]
    ecdfs_marker_size = 5

    # Save variables
    save_filename_absolute_count = "cumbria_devon_premises_type_counts_julia.pdf"
    save_filename_propns = "cumbria_devon_premises_type_propns_julia.pdf"
    save_filename_popn_ecdfs = "cumbria_devon_popn_ecdfs_julia.pdf"
    save_filename_county_locator_map = "cumbria_devon_locator_map_julia.pdf"

    #----------------------------------------------------------------------
    # Load the associated livestock population data
    #----------------------------------------------------------------------

    # For each county, retain premises that have either cattle or sheep
    # recorded as being present

    # Initialise cell to store livestock popn info for each county
    livestock_data_tuple = Array{Array{Int64,2},1}(undef,n_counties)

    # Initialise summary statistic storage vectors
    # Row per county.
    # Col 1: Cattle only; Col 2: Sheep only: Col 3: Both cattle and sheep.
    prem_count_livestock_types = zeros(Int64,3,n_counties)
    propn_prem_per_livestock_composition = zeros(Float64,n_counties,3)

    # Iterate over each county and populate livestock_data_tuple
    for county_itr = 1:n_counties

        # Get the county ID, based on the loop iteration index
        county_ID = county_IDs[county_itr]

        # Get the filename containing the relevant livestock popn data for the
        # county
        if county_ID == 0
            county_livestock_data_filename = "../../../data/dummy_data/livestock_dummy_data.txt"
        elseif county_ID < 10
            county_livestock_data_filename = "../../../data/dummy_data/cumbria_synthetic_data.txt"
            # county_livestock_data_filename = "../../../ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID0$county_ID.txt"
        else
            county_livestock_data_filename = "../../../data/dummy_data/devon_synthetic_data.txt"
            # county_livestock_data_filename = "../../../ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID$county_ID.txt"
        end

        # Read in the data and retain the premises (rows) that have positive
        # counts for cattle and/or sheep. Otherwise, omit those entries.
        prem_livestock_data_all = readdlm(county_livestock_data_filename)
            #     Column breakdown
            #     Cols 1-4: Survey year, CPH  Easting, Northing,
            #     Cols 5-10: Cattle, Sheep, Horses, Pigs, Geese, Chickens,
            #     Cols 11-15 Ducks, Goats, Turkeys, Other Poultry, Deer,
            #     Cols 16-20 Other livestock, Alpacas, Llamas, Bee hives
        cattle_and_sheep_column_idxs = [5,6]
        prem_livestock_data_cattle_and_sheep_only = prem_livestock_data_all[:,cattle_and_sheep_column_idxs]

        # # Check premises lies within specified bounding box
        # if county_ID == 10
        #     # Devon
        #
        #     # Bounding box
        #     bounding_box_left = 212000
        #     bounding_box_right = 340000
        #     bounding_box_bottom = 34960
        #     bounding_box_top = 151210
        #
        #     # Check validity of locations.
        #     # Remove locations outside of bounding box
        #     prem_loc_raw = prem_livestock_data_all[:,[3,4]]
        #     valid_prem_loc_flag = (prem_loc_raw[:,1] .>= bounding_box_left) .&
        #                         (prem_loc_raw[:,1] .<= bounding_box_right) .&
        #                         (prem_loc_raw[:,2] .>= bounding_box_bottom) .&
        #                         (prem_loc_raw[:,2] .<= bounding_box_top)
        # else
        #
        # end

        # Retain those premises that have cattle and/or sheep present AND
        # are within the bounding box
        total_cattle_sheep_on_prem = sum(prem_livestock_data_cattle_and_sheep_only,dims=2)[:]
            # Use [:] to put into vector
        valid_prem_loc_flag = ones(Bool,size(prem_livestock_data_all,1))
        retain_prem_flag = (total_cattle_sheep_on_prem.>0) .& (valid_prem_loc_flag .== true)

        # Assign retained premises location and livestock records to storage cell
        location_cattle_and_sheep_column_idxs = [3;4;5;6]
        livestock_data_tuple[county_itr] = prem_livestock_data_all[retain_prem_flag,location_cattle_and_sheep_column_idxs]

        # Populate livestock type array (absolute counts)
        # Column assignment:
        # Cattle only present, Sheep only present, Both cattle and sheep present.
        livestock_data_prem_with_cattle_or_sheep = prem_livestock_data_cattle_and_sheep_only[retain_prem_flag,:]
        prem_count_livestock_types[1,county_itr] = sum((livestock_data_prem_with_cattle_or_sheep[:,1] .> 0) .& (livestock_data_prem_with_cattle_or_sheep[:,2] .== 0))
        prem_count_livestock_types[2,county_itr] = sum((livestock_data_prem_with_cattle_or_sheep[:,1] .== 0) .& (livestock_data_prem_with_cattle_or_sheep[:,2] .> 0))
        prem_count_livestock_types[3,county_itr] = sum((livestock_data_prem_with_cattle_or_sheep[:,1] .> 0) .& (livestock_data_prem_with_cattle_or_sheep[:,2] .> 0))

        # Error check
        if ( sum(prem_count_livestock_types[:,county_itr]) != sum(retain_prem_flag) )
            error("Sum of premises in prem_count_livestock_types is incorrect for county_itr $county_itr")
        end

        # Populate livestock type array (proportions)
        propn_prem_per_livestock_composition[county_itr,:] = prem_count_livestock_types[:,county_itr]./sum(prem_count_livestock_types[:,county_itr])
    end

    #----------------------------------------------------------------------
    # Absolute counts (cattle only, sheep only, both)
    #----------------------------------------------------------------------

    # Generate bar plot
    plt_bar_counts = groupedbar(absolute_counts_groups, # Group labels on x-axis
                                prem_count_livestock_types, # y-axis
                                linewidth = 1.5,
                                bar_width = 0.8,
                                color = absolute_counts_colourmap,
                                label = absolute_counts_leg_labels,
                                framestyle = :box,
                                bar_position = :dodge,
                                grid = false)

    # # Add values to top of each bar
    x_pos = [0.30 0.7;
             1.30 1.7;
             2.30 2.7]
    for county_itr = 1:n_counties
        for grp_itr = 1:3
            annotate!(x_pos[grp_itr,county_itr], # x position
                      prem_count_livestock_types[grp_itr,county_itr]+100, # y position
                      text(string(prem_count_livestock_types[grp_itr,county_itr]), :black, :centre, plot_fontsize_default)) # Text to add to plot
        end
    end


    # Add legend
    plot!(legend = true,
            legendfontsize = plot_fontsize_default)

    # Set axes labels
    plot!(ylabel = "Number of premises")

    # Set axis limits
    plot!(ylims = (0,2500))

    # Set ticklabels
    plot!(yticks = (0:500:2500, 0:500:2500))

    # Set fontsizes
    plot!(xtickfontsize = plot_fontsize_default,
            ytickfontsize = plot_fontsize_default,
            guidefontsize = plot_fontsize_default)

    # If appropriate, specify the overall plot size
    plot!(size=(fig_width, fig_height))

    # If applicable, save the figure
    if save_fig_flag == true
        savefig(plt_bar_counts, save_filename_absolute_count)
    end

    #----------------------------------------------------------------------
    # Proportions (cattle only, sheep only, both) as stacked bar plot
    #----------------------------------------------------------------------

    # Set group x-axis positions
    X_propns = 1:1:length(propns_groups)

    # Generate stacked bar plot
    plt_bar_propns = groupedbar(X_propns, # Group labels on x-axis
                                propn_prem_per_livestock_composition, # y-axis
                                linewidth = 1.5,
                                bar_width = 0.67,
                                color = propns_colourmap,
                                label = ["Cattle only" "Sheep only" "Cattle and sheep"],
                                framestyle = :box,
                                bar_position = :stack,
                                grid = false)

    # Add legend
    plot!(legend = :outerright,
            legendfontsize = plot_fontsize_default)

    # Set y-axis properties
    plot!(ylabel = "Proportion of premises",
           yticks = (0:0.1:1, ["0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"]),
           yguidefontsize = plot_fontsize_default,
           ytickfontsize = plot_fontsize_default,
           yformatter = :plain)

    # Set x-axis properties
    plot!(xlims = (0.5,length(propns_groups)+0.5),
            xticks = (X_propns, propns_groups),
            xguidefontsize = plot_fontsize_default,
            xtickfontsize = plot_fontsize_default
            )

    # If applicable, save the figure
    if save_fig_flag == true
        savefig(plt_bar_propns, save_filename_propns)
    end

    #----------------------------------------------------------------------
    # Empirical cumulative distribution function (ecdf) of herd and flock sizes
    #----------------------------------------------------------------------

    # Initialise plot
    plt_ecdfs = StatsPlots.plot()

    # Set percentile summary stats for cow herd/sheep flock size to be calculated
    livestock_popn_prctile_vals = [0.5 0.25 0.75 0.975]

    # Extract the cattle herd sizes and sheep flock sizes for each county
    # Compute the ecdfs and add to plot
    for county_itr = 1:n_counties

        # Get IDs of premises where the given livestock type is present
        cattle_herd_present_ID = livestock_data_tuple[county_itr][:,3] .> 0
        sheep_flock_present_ID = livestock_data_tuple[county_itr][:,4] .> 0

        # Construct population vectors limited to premises where the given livestock type was present
        cattle_herd_size_present_vec = livestock_data_tuple[county_itr][cattle_herd_present_ID,3]
        sheep_flock_size_present_vec = livestock_data_tuple[county_itr][sheep_flock_present_ID,4]

        # Print median, IQR and 97.5th percentile information for herd/flock size
        # Only include premises where that flock type was present
        cattle_herd_sum_stats = quantile(cattle_herd_size_present_vec,livestock_popn_prctile_vals)
        sheep_flock_sum_stats = quantile(sheep_flock_size_present_vec,livestock_popn_prctile_vals)
        if county_itr == 1 # Cumbria
            println("Cumbria cattle herd summary stats. Median: $(cattle_herd_sum_stats[1]); IQR: ($(cattle_herd_sum_stats[2]),$(cattle_herd_sum_stats[3])); 97.5th percentile: $(cattle_herd_sum_stats[4])")
            println("Cumbria sheep flock summary stats. Median: $(sheep_flock_sum_stats[1]); IQR: ($(sheep_flock_sum_stats[2]),$(sheep_flock_sum_stats[3])); 97.5th percentile: $(sheep_flock_sum_stats[4])")
        elseif county_itr == 2 # Devon
            println("Devon cattle herd summary stats. Median: $(cattle_herd_sum_stats[1]); IQR: ($(cattle_herd_sum_stats[2]),$(cattle_herd_sum_stats[3])); 97.5th percentile: $(cattle_herd_sum_stats[4])")
            println("Devon sheep flock summary stats. Median: $(sheep_flock_sum_stats[1]); IQR: ($(sheep_flock_sum_stats[2]),$(sheep_flock_sum_stats[3])); 97.5th percentile: $(sheep_flock_sum_stats[4])")
        else
            error("Invalid county_itr value of $county_itr. Expected value of 1 or 2.")
        end

        # Get the ecdfs and x-axis positions
        f_cattle_plot = ecdf(cattle_herd_size_present_vec)
        x_cattle_plot = sort(cattle_herd_size_present_vec)

        f_sheep_plot = ecdf(sheep_flock_size_present_vec)
        x_sheep_plot = sort(sheep_flock_size_present_vec)

        # Add ecdf for cattle herd size to plot
        plot!(x_cattle_plot,y -> f_cattle_plot(y),
                xaxis=:log,
                linestyle = :solid,
                linewidth = 1.5,
                color = ecdfs_colourmap[county_itr],
                markersize = ecdfs_marker_size,
                markerstrokewidth = 1.5,
                markerstrokecolor = ecdfs_colourmap[county_itr],
                framestyle = :box,
                seriestype = :samplemarkers,
                step = 200,
                shape = ecdfs_marker_type[county_itr],
                label = "")

        # Add ecdf for sheep flock size to plot
        plot!(x_sheep_plot,y -> f_sheep_plot(y),
                xaxis=:log,
                linestyle = :dot,
                linewidth = 2.5,
                color = ecdfs_colourmap[county_itr],
                markersize = ecdfs_marker_size,
                markerstrokewidth = 1.5,
                markerstrokecolor = ecdfs_colourmap[county_itr],
                framestyle = :box,
                seriestype = :samplemarkers,
                step = 200,
                shape = ecdfs_marker_type[county_itr],
                label = "")
    end

    # Set labels for legend.
    plot!([0:0.1:0.1],[0:0.1:0.1],
                        linestyle = :solid,
                        linewidth = 1.5,
                        color = ecdfs_colourmap[1],
                        markershape = ecdfs_marker_type[1],
                        markerstrokecolor = ecdfs_colourmap[1],
                        markerstrokewidth = 1.5,
                        markersize = ecdfs_marker_size,
                        label = ecdf_leg_labels[1])
    plot!([0:0.1:0.1],[0:0.1:0.1],
                        linestyle = :dot,
                        linewidth = 2.5,
                        color = ecdfs_colourmap[1],
                        markershape = ecdfs_marker_type[1],
                        markerstrokecolor = ecdfs_colourmap[1],
                        markerstrokewidth = 1.5,
                        markersize = ecdfs_marker_size,
                        label = ecdf_leg_labels[2])
    plot!([0:0.1:0.1],[0:0.1:0.1],
                        linestyle = :solid,
                        linewidth = 1.5,
                        color = ecdfs_colourmap[2],
                        markershape = ecdfs_marker_type[2],
                        markerstrokecolor = ecdfs_colourmap[2],
                        markerstrokewidth = 1.5,
                        markersize = ecdfs_marker_size,
                        label = ecdf_leg_labels[3])
    plot!([0:0.1:0.1],[0:0.1:0.1],
                        linestyle = :dot,
                        linewidth = 2.5,
                        color = ecdfs_colourmap[2],
                        markershape = ecdfs_marker_type[2],
                        markerstrokecolor = ecdfs_colourmap[2],
                        markerstrokewidth = 1.5,
                        markersize = ecdfs_marker_size,
                        label = ecdf_leg_labels[4])
    plot!(legend = :bottomright,
                    legendfontsize = plot_fontsize_ecdfs_leg,
                    xminorticks = true) # Also adds minorticks to log x-axis

    # Set y-axis properties
    plot!(ylabel = "Empirical CDF",
           ylims = (0,1),
           yticks = (0:0.1:1, ["0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"]),
           yguidefontsize = plot_fontsize_ecdfs,
           ytickfontsize = plot_fontsize_ecdfs,
           yformatter = :plain)

    # Set x-axis properties
    plot!(xlabel = "Livestock population",
            xlims = (1,10000),
            xguidefontsize = plot_fontsize_ecdfs,
            xtickfontsize = plot_fontsize_ecdfs
            )

    # Specify the overall plot size
    plot!(size=(fig_width_ecdfs, fig_height_ecdfs))

    # If applicable, save the figure
    if save_fig_flag == true
        savefig(plt_ecdfs, save_filename_popn_ecdfs)
    end


    #---------------------------------------------------------------------------
    # Locator map
    #---------------------------------------------------------------------------

    # Load shapefiles data
    shapefile_table =  Shapefile.Table("../../../data/Shapefiles/Counties_and_Unitary_Authorities_(December_2020)_UK_BGC/Counties_and_Unitary_Authorities_(December_2020)_UK_BGC.shp")
    county_shapes = Shapefile.shapes(shapefile_table)

    # Generate all county outlines, no filled polygons
    plt_map = plot(county_shapes[1:151],
                    color = :white,
                    linecolor = :blue,
                    linewidth = 0.5,
                    grid = false,
                    showaxis = false,
                    xticks = nothing,
                    yticks = nothing,
                    xlims = (100000,700000),
                    label = "",
                    size = (450,500))

    # Shade Cumbria
    plot!(county_shapes[128],
            color = RGB(0.2,0.2,0.2),
            linecolor = :blue,
            linewidth = 0,
            label = "Cumbria")

    # Shade Devon
    plot!(county_shapes[130],
            color = RGB(0.8,0.8,0.8),
            linecolor = :blue,
            linewidth = 0,
            label = "Devon")

    # Add legend
    plot!(legend = :topright,
            legendfontsize = plot_fontsize_map)

    # If applicable, save the figure
    if save_fig_flag == true
        savefig(plt_map, save_filename_county_locator_map)
    end

    return nothing
end
