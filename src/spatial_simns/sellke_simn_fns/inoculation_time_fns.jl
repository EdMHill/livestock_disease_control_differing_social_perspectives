#=
Purpose:
File to house functions to set delay in vaccination becoming effective

Function list:
zero_days_inoculation_fn - No delay in vaccine becoming effective
one_day_inoculation_fn - Single day delay in vaccine becoming effective
uniform_four_to_six_days_inoculation_fn - Equal chance of 4, 5 or 6 day delay.

Date: 3rd November 2021
=#


#= FUNCTION GENERIC OUTLINE

    inoculation_fn(rng,PremNum)

Randomly draw sample from specified discrete probability distribution.
Used to randomly set infection, adherence and testing waiting times.

Inputs: `rng` - random number generator,
        `PremNum` - Number of premises in the simulation
Outputs: `prem_time_to_inoculation` - For each premises, delay in vaccine becoming effective once administered
=#

#-------------------------------------------------------------------------------
### ZERO_DAYS_INOCULATION_FN - NO DELAY IN VACCINE BECOMING EFFECTIVE
#-------------------------------------------------------------------------------
"""
    zero_days_inoculation_fn(rng::AbstractRNG,PremNum::Int64)

Inoculation function: No delay in vaccine becoming effective.

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `PremNum::Int64`: Number of premises.

Outputs:
- `prem_time_to_inoculation::Array{Float64,1}`: Per premises, delay in vaccine becoming effective after being administered.

Location: inoculation\\_time\\_fns.jl
"""
function zero_days_inoculation_fn(rng::AbstractRNG,PremNum::Int64)

    # Return a zero vector
    prem_time_to_inoculation = zeros(Float64,PremNum)

    return prem_time_to_inoculation::Array{Float64,1}
end

#-------------------------------------------------------------------------------
### ONE_DAY_INOCULATION_FN - SINGLE DAY DELAY IN VACCINE BECOMING EFFECTIVE
#-------------------------------------------------------------------------------
"""
    one_day_inoculation_fn(rng::AbstractRNG,PremNum::Int64)

Inoculation function: Single day delay in vaccine becoming effective.

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `PremNum::Int64`: Number of premises.

Outputs:
- `prem_time_to_inoculation::Array{Float64,1}`: Per premises, delay in vaccine becoming effective after being administered.

Location: inoculation\\_time\\_fns.jl
"""
function one_day_inoculation_fn(rng::AbstractRNG,PremNum::Int64)

    # Return a zero vector
    prem_time_to_inoculation = ones(Float64,PremNum)

    return prem_time_to_inoculation::Array{Float64,1}
end


#-------------------------------------------------------------------------------
### UNIFORM_FOUR_TO_SIX_DAYS_INOCULATION_FN - EQUAL CHANCE OF 4, 5 OR 6 DAY DELAY.
#-------------------------------------------------------------------------------
"""
    uniform_four_to_six_days_inoculation_fn(rng::AbstractRNG,PremNum::Int64)

Inoculation function: Equal chance of a 4, 5 or 6 day delay in vaccine becoming effective.

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `PremNum::Int64`: Number of premises.

Outputs:
- `prem_time_to_inoculation::Array{Float64,1}`: Per premises, delay in vaccine becoming effective after being administered.

Location: inoculation\\_time\\_fns.jl
"""
function uniform_four_to_six_days_inoculation_fn(rng::AbstractRNG,PremNum::Int64)

    # Set possible delay times and weights
    inoculation_time_vals = [4.,5.,6.]
    inoculation_time_weights = [1/3,1/3,1/3]

    # Sample the inoculation time for each premises
    prem_time_to_inoculation = sample(rng,inoculation_time_vals,Weights(inoculation_time_weights),PremNum)

    return prem_time_to_inoculation::Array{Float64,1}
end
