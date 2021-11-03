#=
Purpose:
Construct array corresponding to control strategy variable combinations

Julia version: 1.6.3
=#

#-------------------------------------------------------------------------------
# LOAD REQUIRED ENVIRONMENT
#-------------------------------------------------------------------------------

#Set path to directory this file resides in
cd(dirname(@__FILE__))

# Activate the environment
using Pkg
Pkg.activate("../../../../")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
using DelimitedFiles

#-------------------------------------------------------------------------------
# COMPUTE NUMBER OF SCENARIOS
#-------------------------------------------------------------------------------

# Possible reactionary strategies
# Respond if notified infection within specified zone (distance in metres)
reactionary_strat_vals = [-1;collect(1000:1000:10000)]
n_reactionary_strats = length(reactionary_strat_vals)

# Stratification of group proportions
permitted_grp_propn_vals = collect(0:0.05:1)
n_permitted_grp_propn_vals = length(permitted_grp_propn_vals)

# Get total number of scenarios
n_scens = sum(collect(1:n_permitted_grp_propn_vals))*n_reactionary_strats

#-------------------------------------------------------------------------------
# INITIALISE ARRAY
#-------------------------------------------------------------------------------
control_param_val_combination_array = zeros(Float64,n_scens,3)

#-------------------------------------------------------------------------------
# POPULATE ARRAY
#-------------------------------------------------------------------------------
row_assign_idx = 1
for initial_adopter_grp_itr = 1:n_permitted_grp_propn_vals
    # For current iteration of loop, assign propn in initial adopter group to variable
    initial_adopter_grp_val = permitted_grp_propn_vals[initial_adopter_grp_itr]
    for reactionary_grp_itr = 1:((n_permitted_grp_propn_vals-initial_adopter_grp_itr)+1)
        # Index goes up to "(n_permitted_group_propn_vals-initial_adopter_grp_itr)+1" to ensure
        # the sum across the initial adopter and reactionary groups does not exceed 1

        # For current iteration of loop, assign propn in reactionary group to variable
        reactionary_grp_val = permitted_grp_propn_vals[reactionary_grp_itr]

        # Error check
        # If sum across the initial adopter and reactionary groups exceeds 1, exit the programme
        if (initial_adopter_grp_val + reactionary_grp_val) > 1
            error("initial_adopter_grp_val + reactionary_grp_val: $(initial_adopter_grp_val + reactionary_grp_val). Group sum exceeds 1, invalid.")
        end

        for reactionary_strat_itr = 1:n_reactionary_strats
            # # Assign values to array
            control_param_val_combination_array[row_assign_idx,1] = initial_adopter_grp_val
            control_param_val_combination_array[row_assign_idx,2] = reactionary_grp_val
            control_param_val_combination_array[row_assign_idx,3] = reactionary_strat_vals[reactionary_strat_itr]

            # Increment row assignment value
            global row_assign_idx += 1
        end

    end
end

#-------------------------------------------------------------------------------
# SAVE ARRAY TO FILE
#-------------------------------------------------------------------------------
writedlm("generic_model_control_param_val_array.txt",control_param_val_combination_array)
