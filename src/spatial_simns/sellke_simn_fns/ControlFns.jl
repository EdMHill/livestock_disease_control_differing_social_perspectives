#=
Purpose:
File to house functions to carry out controls

Function list:
- Infected premises cull + ring culling (calculate distances on the fly)
- Infected premises cull + ring culling (takes precalculated distances as an input)
- Infected premises cull + ring vaccination (all livestock types)
- Infected premises cull + ring vaccination of cattle only
- Infected premises cull + vaccination based on notified infection occuring within a specified distance (versions for all livestock vaccinated & vaccination of cattle only)

Supporting fns:
- cull_prem_and_update_time_to_cattle_vacc_effective_fn!
- apply_cattle_only_vacc_fn! (For a single premises, apply cattle vaccination strategy + its effect)
- deploy_cattle_only_vacc_fn! (For a single premises, apply cattle vaccination strategy with time lag until effect)
- activate_cattle_only_vacc_fn! (For a single premises, apply effect of cattle vaccination strategy only)
- cull_prem_and_update_time_to_all_livestock_vacc_effective_fn!
- apply_all_livestock_vacc_fn! (For a single premises, apply livestock vaccination strategy + its effect)
- deploy_all_livestock_vacc_fn! (For a single premises, apply livestock vaccination strategy with time lag until effect)
- activate_all_livestock_vacc_fn! (For a single premises, apply effect of livestock vaccination strategy)

Julia version: 1.6.3
Date: 3rd November 2021
=#

#-------------------------------------------------------------------------------
# MAIN FUNCTIONS
#-------------------------------------------------------------------------------

"""
    CullIPsAndRingCullFn!(PremNum::Int64,
                            PremStatus::Array{Float64,1},
                            PremLoc_xVals::Array{Float64,1},
                            PremLoc_yVals::Array{Float64,1},
                            EpiParamVals::Array{Float64,1},
                            ControlParamVals::Array{Float64,1},
                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                            CoordType::Int64,
                            CullPremDuringCurrentTimestep::Array{Bool,1})

Apply control: Infected premises cull + ring culling (calculate distances on the fly)

Inputs:
- `PremNum::Int64`: Number of premises in landscape.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)].
- `ControlParamVals::Array{Float64,1}`: Extent of control enforced.
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.

Outputs: None \n
Location: ControlFns.jl
"""
function CullIPsAndRingCullFn!(PremNum::Int64,
                            PremStatus::Array{Float64,1},
                            PremLoc_xVals::Array{Float64,1},
                            PremLoc_yVals::Array{Float64,1},
                            EpiParamVals::Array{Float64,1},
                            ControlParamVals::Array{Float64,1},
                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                            CoordType::Int64,
                            CullPremDuringCurrentTimestep::Array{Bool,1})

    #Disaggregate relevant variables from EpiParamVals
    RemovalTime = EpiParamVals[3]::Float64

    #Disaggregate relevant variables from ControlParamVals
    RingCull = ControlParamVals[1]::Float64

    #Find and update premises to be culled
    m=findall(PremStatus.>RemovalTime)
    for ii = 1:length(m)
        #Get ID of premises to be culled
        PremToCullID = m[ii]

        #Cull IP
        PremStatus[PremToCullID] = -1.

        #Update premises culled during current timestep vector
        CullPremDuringCurrentTimestep[PremToCullID] = 1

        #Ascertain locations to undergo ring cull
        if RingCull > 0
            D = zeros(PremNum)
            if CoordType == 1 #Cartesian (metres unit)
                for DistanceCheckItr = 1:PremNum
                    xDiff = PremLoc_xVals[PremToCullID] - PremLoc_xVals[DistanceCheckItr]
                    yDiff = PremLoc_yVals[PremToCullID] - PremLoc_yVals[DistanceCheckItr]
                    D[DistanceCheckItr] = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                end
            elseif CoordType == 2 #Cartesian (km unit)
                for DistanceCheckItr = 1:PremNum
                    xDiff = PremLoc_xVals[PremToCullID] - PremLoc_xVals[DistanceCheckItr]
                    yDiff = PremLoc_yVals[PremToCullID] - PremLoc_yVals[DistanceCheckItr]
                    D[DistanceCheckItr] = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                end
            elseif CoordType == 3 #LatLong
                for DistanceCheckItr = 1:PremNum
                    D[DistanceCheckItr] = GreatCircleDistance(PremLoc_yVals[PremToCullID], PremLoc_xVals[PremToCullID],
                                        PremLoc_yVals[DistanceCheckItr], PremLoc_xVals[DistanceCheckItr])
                end
            end

            #Find premises within ring cull region
            #Perform ring cull
            for D_Itr = 1:PremNum
                if (PremStatus[D_Itr] != -1) #Check premises not already culled
                    if (D[D_Itr]<=RingCull) #Check if premises is within control radius
                        PremStatus[D_Itr] = -1.
                        CullPremDuringCurrentTimestep[D_Itr] = 1
                    end
                end
            end
            #PremStatus[D.<RingCull] .= -1
        end
    end
    return nothing
    # return PremStatus::Array{Float64,1},
    #         CullPremDuringCurrentTimestep::Array{Bool,1}
end

"""
    CullIPsAndRingCullFn_WithPreCalcDist!(PremNum::Int64,
                                            PremStatus::Array{Float64,1},
                                            PremLoc_xVals::Array{Float64,1},
                                            PremLoc_yVals::Array{Float64,1},
                                            EpiParamVals::Array{Float64,1},
                                            ControlParamVals::Array{Float64,1},
                                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                            CoordType::Int64,
                                            CullPremDuringCurrentTimestep::Array{Bool,1})

Apply control: Infected premises cull + ring culling (takes precalculated distances as an input)

Inputs:
- `PremNum::Int64`: Number of premises in landscape.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)].
- `ControlParamVals::Array{Float64,1}`: Extent of control enforced.
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.

Outputs: None \n
Location: ControlFns.jl
"""
function CullIPsAndRingCullFn_WithPreCalcDist!(PremNum::Int64,
                                                PremStatus::Array{Float64,1},
                                                PremLoc_xVals::Array{Float64,1},
                                                PremLoc_yVals::Array{Float64,1},
                                                EpiParamVals::Array{Float64,1},
                                                ControlParamVals::Array{Float64,1},
                                                LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                                CoordType::Int64,
                                                CullPremDuringCurrentTimestep::Array{Bool,1})

    #Disaggregate relevant variables from EpiParamVals
    RemovalTime = EpiParamVals[3]::Float64

    #Disaggregate relevant variables from ControlParamVals
    RingCull = ControlParamVals[1]::Float64

    #Find and update premises to be culled
    m=findall(PremStatus.>RemovalTime)
    for ii = 1:length(m)
        #Get ID of premises to be culled
        PremToCullID = m[ii]

        #Cull IP
        PremStatus[PremToCullID] = -1.

        #Update premises culled during current timestep vector
        CullPremDuringCurrentTimestep[PremToCullID] = 1

        #Cull premises within ring cull region
        if RingCull > 0
            PremStatus[LinkedPremToControlIdxs[PremToCullID]] .= -1.
        end
    end
    return nothing
    # return PremStatus::Array{Float64,1},
    #         CullPremDuringCurrentTimestep::Array{Bool,1}
end

"""
    CullIPsAndRingVaccFn!(PremNum::Int64,
                            PremStatus::Array{Float64,1},
                            PremVaccStatus::Array{Int64,1},
                            SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                            PremSuscept::Array{Float64},
                            PremTransmiss::Array{Float64},
                            PremLoc_xVals::Array{Float64,1},
                            PremLoc_yVals::Array{Float64,1},
                            EpiParamVals::Array{Float64,1},
                            ControlParamVals::Array{Float64,1},
                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                            CoordType::Int64,
                            CullPremDuringCurrentTimestep::Array{Bool,1},
                            VaccPremDuringCurrentTimestep::Array{Bool,1},
                            VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1})

Apply control: Infected premises cull + ring vaccination (all livestock types, though not vaccinated if reporting infection).

Inputs:
- `PremNum::Int64`: Number of premises in landscape.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)].
- `ControlParamVals::Array{Float64,1}`: Extent of control enforced.
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.

Outputs: None \n
Location: ControlFns.jl
"""
function CullIPsAndRingVaccFn!(PremNum::Int64,
                            PremStatus::Array{Float64,1},
                            PremVaccStatus::Array{Int64,1},
                            SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                            PremSuscept::Array{Float64},
                            PremTransmiss::Array{Float64},
                            PremLoc_xVals::Array{Float64,1},
                            PremLoc_yVals::Array{Float64,1},
                            EpiParamVals::Array{Float64,1},
                            ControlParamVals::Array{Float64,1},
                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                            CoordType::Int64,
                            CullPremDuringCurrentTimestep::Array{Bool,1},
                            VaccPremDuringCurrentTimestep::Array{Bool,1},
                            VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1})

    #Disaggregate relevant variables from EpiParamVals
    DetectionTime = EpiParamVals[2]::Float64
    RemovalTime = EpiParamVals[3]::Float64

    #Disaggregate relevant variables from ControlParamVals
    VaccEff = ControlParamVals[1]::Float64
    RingVaccRadius = ControlParamVals[2]::Float64

    #Find and update premises to be culled
    m=findall(PremStatus.>RemovalTime)
    for ii = 1:length(m)
        #Get ID of premises to be culled
        PremToCullID = m[ii]

        #Cull IP
        PremStatus[PremToCullID] = -1.

        #Update premises culled during current timestep vector
        CullPremDuringCurrentTimestep[PremToCullID] = 1

        #Ascertain locations to undergo ring cull
        if RingVaccRadius > 0
            D = zeros(PremNum)
            if CoordType == 1 #Cartesian (metres unit)
                for DistanceCheckItr = 1:PremNum
                    xDiff = PremLoc_xVals[PremToCullID] - PremLoc_xVals[DistanceCheckItr]
                    yDiff = PremLoc_yVals[PremToCullID] - PremLoc_yVals[DistanceCheckItr]
                    D[DistanceCheckItr] = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                end
            elseif CoordType == 2 #Cartesian (km unit)
                for DistanceCheckItr = 1:PremNum
                    xDiff = PremLoc_xVals[PremToCullID] - PremLoc_xVals[DistanceCheckItr]
                    yDiff = PremLoc_yVals[PremToCullID] - PremLoc_yVals[DistanceCheckItr]
                    D[DistanceCheckItr] = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                end
            elseif CoordType == 3 #LatLong
                for DistanceCheckItr = 1:PremNum
                    D[DistanceCheckItr] = GreatCircleDistance(PremLoc_yVals[PremToCullID], PremLoc_xVals[PremToCullID],
                                        PremLoc_yVals[DistanceCheckItr], PremLoc_xVals[DistanceCheckItr])
                end
            end

            #Find premises within ring cull region
            #Perform ring vaccination if not previously vaccinated AND premises not already culled AND not reported infection
            #Vaccine ineffective if already infected
            for D_Itr = 1:PremNum
                # if (D[D_Itr]<=RingVaccRadius)
                #     println("Disease status for vaccine deployed: $(PremStatus[D_Itr])")
                #     println("LogicCheck: $(0<=PremStatus[D_Itr]<DetectionTime)")
                #     println("LogicCheck: $((0<=PremStatus[D_Itr]) && (PremStatus[D_Itr]<DetectionTime))")
                #     println("PremVaccStatus: $(PremVaccStatus[D_Itr])")
                # end
                if (D[D_Itr]<=RingVaccRadius) && (0<=PremStatus[D_Itr]<DetectionTime) && (PremVaccStatus[D_Itr] == 0)

                    #Vaccinate premises, update PremVaccStatus & SpeciesGroupVaccStatusByPrem
                    PremVaccStatus[D_Itr] = 1
                    SpeciesGroupVaccStatusByPrem[D_Itr,:] .= 1 #All species are vaccinated

                    #Update premises vaccinated during current timestep vector
                    VaccPremDuringCurrentTimestep[D_Itr] = 1

                    #Update vaccinated premises susceptibility and transmissibility
                    #Vaccine only successful in modifying susceptiblity & transmissiblity if susceptible
                    if PremStatus[D_Itr] == 0
                        PremSuscept[D_Itr] = PremSuscept[D_Itr].*(1-VaccEff)
                        PremTransmiss[D_Itr] = PremTransmiss[D_Itr].*(1-VaccEff)
                    end
                end
            end
        end
    end

    return nothing
end

"""
    CullIPsAndRingVaccCattleOnlyFn!(PremNum::Int64,
                                    PremStatus::Array{Float64,1},
                                    PremVaccStatus::Array{Int64,1},
                                    SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                    PremSuscept::Array{Float64},
                                    PremTransmiss::Array{Float64},
                                    PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                    PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                    PremLoc_xVals::Array{Float64,1},
                                    PremLoc_yVals::Array{Float64,1},
                                    EpiParamVals::Array{Float64,1},
                                    ControlParamVals::Array{Float64,1},
                                    LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                    CoordType::Int64,
                                    CullPremDuringCurrentTimestep::Array{Bool,1},
                                    VaccPremDuringCurrentTimestep::Array{Bool,1},
                                    VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                    delta_t::Float64)

Apply control: Infected premises cull + ring vaccination of cattle only (not vaccinated if reporting infection).

Inputs:
- `PremNum::Int64`: Number of premises in landscape.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)].
- `ControlParamVals::Array{Float64,1}`: Extent of control enforced.
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.
- `delta_t::Float64`: Timestep increment.

Outputs: None \n
Location: ControlFns.jl
"""
function CullIPsAndRingVaccCattleOnlyFn!(PremNum::Int64,
                                            PremStatus::Array{Float64,1},
                                            PremVaccStatus::Array{Int64,1},
                                            SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                            PremSuscept::Array{Float64},
                                            PremTransmiss::Array{Float64},
                                            PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                            PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                            PremLoc_xVals::Array{Float64,1},
                                            PremLoc_yVals::Array{Float64,1},
                                            EpiParamVals::Array{Float64,1},
                                            ControlParamVals::Array{Float64,1},
                                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                            CoordType::Int64,
                                            CullPremDuringCurrentTimestep::Array{Bool,1},
                                            VaccPremDuringCurrentTimestep::Array{Bool,1},
                                            VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                            delta_t::Float64)

    #Disaggregate relevant variables from EpiParamVals
    DetectionTime = EpiParamVals[2]::Float64
    RemovalTime = EpiParamVals[3]::Float64

    #Disaggregate relevant variables from ControlParamVals
    VaccEff = ControlParamVals[1]::Float64
    RingVaccRadius = ControlParamVals[2]::Float64

    #Find and update premises to be culled
    m=findall(PremStatus.>RemovalTime)
    for ii = 1:length(m)
        #Get ID of premises to be culled
        PremToCullID = m[ii]

        #Cull IP
        PremStatus[PremToCullID] = -1.

        #Update premises culled during current timestep vector
        CullPremDuringCurrentTimestep[PremToCullID] = 1

        #Ascertain locations to undergo ring cull
        if RingVaccRadius > 0
            D = zeros(PremNum)
            if CoordType == 1 #Cartesian (metres unit)
                for DistanceCheckItr = 1:PremNum
                    xDiff = PremLoc_xVals[PremToCullID] - PremLoc_xVals[DistanceCheckItr]
                    yDiff = PremLoc_yVals[PremToCullID] - PremLoc_yVals[DistanceCheckItr]
                    D[DistanceCheckItr] = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                end
            elseif CoordType == 2 #Cartesian (km unit)
                for DistanceCheckItr = 1:PremNum
                    xDiff = PremLoc_xVals[PremToCullID] - PremLoc_xVals[DistanceCheckItr]
                    yDiff = PremLoc_yVals[PremToCullID] - PremLoc_yVals[DistanceCheckItr]
                    D[DistanceCheckItr] = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                end
            elseif CoordType == 3 #LatLong
                for DistanceCheckItr = 1:PremNum
                    D[DistanceCheckItr] = GreatCircleDistance(PremLoc_yVals[PremToCullID], PremLoc_xVals[PremToCullID],
                                        PremLoc_yVals[DistanceCheckItr], PremLoc_xVals[DistanceCheckItr])
                end
            end

            #Find premises within ring cull region
            #Perform ring vaccination if not previously vaccinated AND premises not already culled AND not reported infection
            #Vaccine ineffective if already infected
            for D_Itr = 1:PremNum
                # if (D[D_Itr]<=RingVaccRadius)
                #     println("Disease status for vaccine deployed: $(PremStatus[D_Itr])")
                #     println("LogicCheck: $(0<=PremStatus[D_Itr]<DetectionTime)")
                #     println("LogicCheck: $((0<=PremStatus[D_Itr]) && (PremStatus[D_Itr]<DetectionTime))")
                #     println("PremVaccStatus: $(PremVaccStatus[D_Itr])")
                # end
                if (D[D_Itr]<=RingVaccRadius) && (0<=PremStatus[D_Itr]<DetectionTime) && (PremVaccStatus[D_Itr] == 0)

                    #Vaccinate premises, update PremVaccStatus & SpeciesGroupVaccStatusByPrem
                    PremVaccStatus[D_Itr] = 1
                    SpeciesGroupVaccStatusByPrem[D_Itr,1] = 1 #Only cattle are vaccinated (first column of livestock type)

                    #Update premises vaccinated during current timestep vector
                    VaccPremDuringCurrentTimestep[D_Itr] = 1

                    #Update vaccinated premises susceptibility and transmissibility for cattle. Modify premises-level values.
                    #Vaccine only successful in modifying susceptiblity & transmissiblity if susceptible
                    if PremStatus[D_Itr] == 0
                        PremSuscept_SpeciesBreakdown[D_Itr,1] = PremSuscept_SpeciesBreakdown[D_Itr,1].*(1-VaccEff)
                        PremTransmiss_SpeciesBreakdown[D_Itr,1] = PremTransmiss_SpeciesBreakdown[D_Itr,1].*(1-VaccEff)
                        PremSuscept[D_Itr] = sum(PremSuscept_SpeciesBreakdown[D_Itr,:])
                        PremTransmiss[D_Itr] = sum(PremTransmiss_SpeciesBreakdown[D_Itr,:])
                    end
                end
            end
        end
    end

    return nothing
    # return PremStatus::Array{Float64,1},
    #         PremVaccStatus::Array{Int64,1},
    #         SpeciesGroupVaccStatusByPrem::Array{Int64,2},
    #         PremSuscept::Array{Float64},
    #         PremTransmiss::Array{Float64},
    #         CullPremDuringCurrentTimestep::Array{Bool,1},
    #         VaccPremDuringCurrentTimestep::Array{Bool,1}
end

"""
    CullIPsAndDistanceBasedVaccFn!(PremNum::Int64,
                                    PremStatus::Array{Float64,1},
                                    PremVaccStatus::Array{Int64,1},
                                    SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                    PremSuscept::Array{Float64},
                                    PremTransmiss::Array{Float64},
                                    PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                    PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                    PremLoc_xVals::Array{Float64,1},
                                    PremLoc_yVals::Array{Float64,1},
                                    EpiParamVals::Array{Float64,1},
                                    ControlParamVals::Array{Any,1},
                                    LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                    CoordType::Int64,
                                    CullPremDuringCurrentTimestep::Array{Bool,1},
                                    VaccPremDuringCurrentTimestep::Array{Bool,1},
                                    VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                    delta_t::Float64)

Apply control: Infected premises cull + vaccination of all livestock types in response to notification of infection within a specified distance.

Inputs:
- `PremNum::Int64`: Number of premises in landscape.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)].
- `ControlParamVals::Array{Float64,1}`: Extent of control enforced.
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.
- `delta_t::Float64`: Timestep increment.

Outputs: None \n
Location: ControlFns.jl
"""
function CullIPsAndDistanceBasedVaccFn!(PremNum::Int64,
                                        PremStatus::Array{Float64,1},
                                        PremVaccStatus::Array{Int64,1},
                                        SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                        PremSuscept::Array{Float64},
                                        PremTransmiss::Array{Float64},
                                        PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                        PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                        PremLoc_xVals::Array{Float64,1},
                                        PremLoc_yVals::Array{Float64,1},
                                        EpiParamVals::Array{Float64,1},
                                        ControlParamVals::Array{Any,1},
                                        LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                        CoordType::Int64,
                                        CullPremDuringCurrentTimestep::Array{Bool,1},
                                        VaccPremDuringCurrentTimestep::Array{Bool,1},
                                        VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                        delta_t::Float64)

    #---------------------------------------------------------------------------
    # DISAGGREGATE RELEVANT VARIABLES
    #---------------------------------------------------------------------------

    #Disaggregate relevant variables from EpiParamVals
    DetectionTime = EpiParamVals[2]::Float64
    RemovalTime = EpiParamVals[3]::Float64

    #Disaggregate vaccination variables from ControlParamVals
    VaccEff = ControlParamVals[1]::Float64

    # Dsiaggregate distance threshold measure
    DistThresVal = ControlParamVals[4]::Float64
        # If infection notified within the given distance, that premises vaccinates

    # Assign vector denoting the phase a premises may undergo vaccination
    # 0 - Not ever vaccinated
    # 1 - Vaccinated prior to simulation start
    # 2 - Vaccinated when risk based measure is satisfied
    PremVaccStage = ControlParamVals[5]::Vector{Int64}

    # Disaggregate time to inoculation variables
    prem_time_to_inoculation = ControlParamVals[7]::Vector{Float64}

    #---------------------------------------------------------------------------
    # FIND AND UPDATE PREMISES TO BE CULLED
    # FIND AND UPDATE TIME TO VACCINE BECOMING EFFECTIVE
    #---------------------------------------------------------------------------
    cull_prem_and_update_time_to_all_livestock_vacc_effective_fn!(PremNum,
                                                                    RemovalTime,
                                                                    PremStatus,
                                                                    PremVaccStatus,
                                                                    CullPremDuringCurrentTimestep,
                                                                    VaccBecomesEffectiveDuringCurrentTimestep,
                                                                    prem_time_to_inoculation,
                                                                    PremSuscept,
                                                                    PremTransmiss,
                                                                    PremSuscept_SpeciesBreakdown,
                                                                    PremTransmiss_SpeciesBreakdown,
                                                                    VaccEff,
                                                                    delta_t)

    #---------------------------------------------------------------------------
    # EVALUATE IF VACCINATION OCCURS
    # FOR HOLDINGS NOT PREVIOUSLY VACCINATED, HAS A NOTIFIED INFECTION OCCURRED ON
    # CURRENT TIMESTEP WITHIN THE THRESHOLD DISTANCE?
    #---------------------------------------------------------------------------

    # NEED PREMISES NOTIFYING INFECTION ON CURRENT TIMESTEP
    newly_notified_infection_prem_IDs = findall(PremStatus.==(DetectionTime+delta_t))
    n_newly_notified_infection_prem_IDs = length(newly_notified_infection_prem_IDs)

    # If there are newly notified infections, enter loop to check whether
    # any premises vaccinate
    if n_newly_notified_infection_prem_IDs > 0

        # Iterate over each premises
        #   - Check if unvaccinated & vaccinates based on risk
        #   - Check distance to premises now reporting infection
        #   - If within the threshold distance, apply vaccination
        for prem_itr = 1:PremNum

            # Check premises is unvaccinated, part of the group that vaccinates in response
            # to risk measure and not had notified infection
            if ( (PremVaccStatus[prem_itr] == 0) &&
                 (PremVaccStage[prem_itr] == 2) &&
                 (0<=PremStatus[prem_itr]<DetectionTime) )

                # Iterate over each premises with newly notified infection
                for newly_notified_prem_itr = 1:n_newly_notified_infection_prem_IDs

                    # Get index of the newly notified premises in this iteration of for loop
                    newly_notified_prem_ID = newly_notified_infection_prem_IDs[newly_notified_prem_itr]

                    # Check distance to premises now reporting infection
                    if CoordType == 1 #Cartesian co-ords (metres)
                        dist_between_prem =  eucl_distance(PremLoc_xVals[newly_notified_prem_ID],
                                            PremLoc_yVals[newly_notified_prem_ID],
                                            PremLoc_xVals[prem_itr],
                                            PremLoc_yVals[prem_itr])
                    elseif CoordType == 2 #Cartesian co-ords (metres)
                        dist_between_prem = eucl_distance_ConvertToMetres(PremLoc_xVals[newly_notified_prem_ID],
                                                            PremLoc_yVals[newly_notified_prem_ID],
                                                            PremLoc_xVals[prem_itr],
                                                            PremLoc_yVals[prem_itr])
                    elseif CoordType == 3 #Lat/Long co-ords
                        println("In LatLong loop!")
                        dist_between_prem = GreatCircleDistance(PremLoc_yVals[newly_notified_prem_ID], PremLoc_xVals[newly_notified_prem_ID],  #lat1, lon1
                                                            PremLoc_yVals[prem_itr], PremLoc_xVals[prem_itr]) #lat2, lon2
                    end

                    # If within the threshold distance, apply vaccination
                    if dist_between_prem <= DistThresVal
                        apply_all_livestock_vacc_fn!(prem_itr,
                                                        PremNum,
                                                        DetectionTime,
                                                        PremStatus,
                                                        PremVaccStatus,
                                                        SpeciesGroupVaccStatusByPrem,
                                                        PremSuscept,
                                                        PremTransmiss,
                                                        PremSuscept_SpeciesBreakdown,
                                                        PremTransmiss_SpeciesBreakdown,
                                                        VaccPremDuringCurrentTimestep,
                                                        VaccBecomesEffectiveDuringCurrentTimestep,
                                                        VaccEff,
                                                        PremVaccStage,
                                                        prem_time_to_inoculation)

                        # Break out the inner for loop.
                        # No need to check other newly notified premises
                        break
                    end
                end
            end
        end
    end

    return nothing
end

"""
    CullIPsAndDistanceBasedCattleVaccFn!(PremNum::Int64,
                                            PremStatus::Array{Float64,1},
                                            PremVaccStatus::Array{Int64,1},
                                            SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                            PremSuscept::Array{Float64},
                                            PremTransmiss::Array{Float64},
                                            PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                            PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                            PremLoc_xVals::Array{Float64,1},
                                            PremLoc_yVals::Array{Float64,1},
                                            EpiParamVals::Array{Float64,1},
                                            ControlParamVals::Array{Any,1},
                                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                            CoordType::Int64,
                                            CullPremDuringCurrentTimestep::Array{Bool,1},
                                            VaccPremDuringCurrentTimestep::Array{Bool,1},
                                            VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                            delta_t::Float64)

Apply control: Infected premises cull + vaccination of cattle in response to notification of infection within a specified distance.

Inputs:
- `PremNum::Int64`: Number of premises in landscape.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)].
- `ControlParamVals::Array{Float64,1}`: Extent of control enforced.
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.
- `delta_t::Float64`: Timestep increment.

Outputs: None \n
Location: ControlFns.jl
"""
function CullIPsAndDistanceBasedCattleVaccFn!(PremNum::Int64,
                                                PremStatus::Array{Float64,1},
                                                PremVaccStatus::Array{Int64,1},
                                                SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                                PremSuscept::Array{Float64},
                                                PremTransmiss::Array{Float64},
                                                PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                                PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                                PremLoc_xVals::Array{Float64,1},
                                                PremLoc_yVals::Array{Float64,1},
                                                EpiParamVals::Array{Float64,1},
                                                ControlParamVals::Array{Any,1},
                                                LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                                CoordType::Int64,
                                                CullPremDuringCurrentTimestep::Array{Bool,1},
                                                VaccPremDuringCurrentTimestep::Array{Bool,1},
                                                VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                                delta_t::Float64)

    #---------------------------------------------------------------------------
    # DISAGGREGATE RELEVANT VARIABLES
    #---------------------------------------------------------------------------

    #Disaggregate relevant variables from EpiParamVals
    DetectionTime = EpiParamVals[2]::Float64
    RemovalTime = EpiParamVals[3]::Float64

    #Disaggregate vaccination variables from ControlParamVals
    VaccEff = ControlParamVals[1]::Float64

    # Dsiaggregate distance threshold measure
    DistThresVal = ControlParamVals[4]::Float64
        # If infection notified within the given distance, that premises vaccinates

    # Assign vector denoting the phase a premises may undergo vaccination
    # 0 - Not ever vaccinated
    # 1 - Vaccinated prior to simulation start
    # 2 - Vaccinated when risk based measure is satisfied
    PremVaccStage = ControlParamVals[5]::Vector{Int64}

    # Disaggregate time to inoculation variables
    prem_time_to_inoculation = ControlParamVals[7]::Vector{Float64}

    #---------------------------------------------------------------------------
    # FIND AND UPDATE PREMISES TO BE CULLED
    # FIND AND UPDATE TIME TO VACCINE BECOMING EFFECTIVE
    #---------------------------------------------------------------------------
    cull_prem_and_update_time_to_cattle_vacc_effective_fn!(PremNum,
                                                            RemovalTime,
                                                            PremStatus,
                                                            PremVaccStatus,
                                                            CullPremDuringCurrentTimestep,
                                                            VaccBecomesEffectiveDuringCurrentTimestep,
                                                            prem_time_to_inoculation,
                                                            PremSuscept,
                                                            PremTransmiss,
                                                            PremSuscept_SpeciesBreakdown,
                                                            PremTransmiss_SpeciesBreakdown,
                                                            VaccEff,
                                                            delta_t)

    #---------------------------------------------------------------------------
    # EVALUATE IF VACCINATION OCCURS
    # FOR HOLDINGS NOT PREVIOUSLY VACCINATED, HAS A NOTIFIED INFECTION OCCURRED ON
    # CURRENT TIMESTEP WITHIN THE THRESHOLD DISTANCE?
    #---------------------------------------------------------------------------

    # NEED PREMISES NOTIFYING INFECTION ON CURRENT TIMESTEP
    newly_notified_infection_prem_IDs = findall(PremStatus.==(DetectionTime+delta_t))
    n_newly_notified_infection_prem_IDs = length(newly_notified_infection_prem_IDs)

    # If there are newly notified infections, enter loop to check whether
    # any premises vaccinate
    if n_newly_notified_infection_prem_IDs > 0

        # Iterate over each premises
        #   - Check if unvaccinated & vaccinates based on risk
        #   - Check distance to premises now reporting infection
        #   - If within the threshold distance, apply vaccination
        for prem_itr = 1:PremNum

            # Check premises is unvaccinated, part of the group that vaccinates in response
            # to risk measure and not had notified infection
            if ( (PremVaccStatus[prem_itr] == 0) &&
                 (PremVaccStage[prem_itr] == 2) &&
                 (0<=PremStatus[prem_itr]<DetectionTime) )

                # Iterate over each premises with newly notified infection
                for newly_notified_prem_itr = 1:n_newly_notified_infection_prem_IDs

                    # Get index of the newly notified premises in this iteration of for loop
                    newly_notified_prem_ID = newly_notified_infection_prem_IDs[newly_notified_prem_itr]

                    # Check distance to premises now reporting infection
                    if CoordType == 1 #Cartesian co-ords (metres)
                        dist_between_prem =  eucl_distance(PremLoc_xVals[newly_notified_prem_ID],
                                            PremLoc_yVals[newly_notified_prem_ID],
                                            PremLoc_xVals[prem_itr],
                                            PremLoc_yVals[prem_itr])
                    elseif CoordType == 2 #Cartesian co-ords (metres)
                        dist_between_prem = eucl_distance_ConvertToMetres(PremLoc_xVals[newly_notified_prem_ID],
                                                            PremLoc_yVals[newly_notified_prem_ID],
                                                            PremLoc_xVals[prem_itr],
                                                            PremLoc_yVals[prem_itr])
                    elseif CoordType == 3 #Lat/Long co-ords
                        println("In LatLong loop!")
                        dist_between_prem = GreatCircleDistance(PremLoc_yVals[newly_notified_prem_ID], PremLoc_xVals[newly_notified_prem_ID],  #lat1, lon1
                                                            PremLoc_yVals[prem_itr], PremLoc_xVals[prem_itr]) #lat2, lon2
                    end

                    # If within the threshold distance, apply vaccination
                    if dist_between_prem <= DistThresVal
                        apply_cattle_only_vacc_fn!(prem_itr,
                                                    PremNum,
                                                    DetectionTime,
                                                    PremStatus,
                                                    PremVaccStatus,
                                                    SpeciesGroupVaccStatusByPrem,
                                                    PremSuscept,
                                                    PremTransmiss,
                                                    PremSuscept_SpeciesBreakdown,
                                                    PremTransmiss_SpeciesBreakdown,
                                                    VaccPremDuringCurrentTimestep,
                                                    VaccBecomesEffectiveDuringCurrentTimestep,
                                                    VaccEff,
                                                    PremVaccStage,
                                                    prem_time_to_inoculation)

                        # Break out the inner for loop.
                        # No need to check other newly notified premises
                        break
                    end
                end
            end
        end
    end

    return nothing
end


#---------------------------------------------------------------------------
# SUPPORTING FUNCTIONS
#---------------------------------------------------------------------------


"""
    cull_prem_and_update_time_to_cattle_vacc_effective_fn!(PremNum::Int64,
                                                            RemovalTime::Float64,
                                                            PremStatus::Array{Float64,1},
                                                            PremVaccStatus::Array{Int64,1},
                                                            CullPremDuringCurrentTimestep::Array{Bool,1},
                                                            VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                                            prem_time_to_inoculation::Array{Float64,1},
                                                            PremSuscept::Array{Float64},
                                                            PremTransmiss::Array{Float64},
                                                            PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                                            PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                                            VaccEff::Float64,
                                                            delta_t::Float64)

Find and update: (i) premises to be culled; (ii) time to vaccine becoming effective (vaccine administered to cattle only).

Inputs:
- `PremNum::Int64`: Number of premises in landscape.
- `RemovalTime::Float64`: Time from infection until livestock on premises are culled.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.
- `prem_time_to_inoculation::Array{Float64,1}`: Time for vaccine to become effective. Entry per premises.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `VaccEff::Float64`: Efficacy of vaccine.
- `delta_t::Float64`: Timestep.

Outputs: None \n
Location: ControlFns.jl
"""
function cull_prem_and_update_time_to_cattle_vacc_effective_fn!(PremNum::Int64,
                                                                RemovalTime::Float64,
                                                                PremStatus::Array{Float64,1},
                                                                PremVaccStatus::Array{Int64,1},
                                                                CullPremDuringCurrentTimestep::Array{Bool,1},
                                                                VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                                                prem_time_to_inoculation::Array{Float64,1},
                                                                PremSuscept::Array{Float64},
                                                                PremTransmiss::Array{Float64},
                                                                PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                                                PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                                                VaccEff::Float64,
                                                                delta_t::Float64
                                                                )
    for prem_idx = 1:PremNum
        # If applicable, find and update premises to be culled
        if PremStatus[prem_idx]>RemovalTime
            #Cull IP
            PremStatus[prem_idx] = -1.

            #Update premises culled during current timestep vector
            CullPremDuringCurrentTimestep[prem_idx] = 1
        end

        # If applicable, find and update time to vaccine becoming effective
        if PremVaccStatus[prem_idx] == 1 # Find premises that are vaccinated

            # Check if time to inoculation is above zero.
            if prem_time_to_inoculation[prem_idx] > 0
                # If satisfied, decrease by timestep
                prem_time_to_inoculation[prem_idx] -= delta_t

                # Error check
                # prem_time_to_inoculation[prem_idx] should not decrease below zero
                if prem_time_to_inoculation[prem_idx] < 0
                    error("prem_time_to_inoculation[$prem_idx]: $(prem_time_to_inoculation[prem_idx]). Invalid.")
                end

                # If time to inoculation reaches zero, apply effect of vaccination (if valid)
                if prem_time_to_inoculation[prem_idx] == 0

                    # Modify premises-level values for susceptibility and transmissibility for cattle.
                    # Vaccine only successful in modifying susceptiblity & transmissiblity if susceptible
                    activate_cattle_only_vacc_fn!(prem_idx,
                                                    PremStatus,
                                                    PremSuscept,
                                                    PremTransmiss,
                                                    PremSuscept_SpeciesBreakdown,
                                                    PremTransmiss_SpeciesBreakdown,
                                                    VaccEff,
                                                    VaccBecomesEffectiveDuringCurrentTimestep)
                end
            end
        end
    end
end

"""
     apply_cattle_only_vacc_fn!(prem_idx::Int64,
                                PremNum::Int64,
                                DetectionTime::Float64,
                                PremStatus::Array{Float64,1},
                                PremVaccStatus::Array{Int64,1},
                                SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                PremSuscept::Array{Float64},
                                PremTransmiss::Array{Float64},
                                PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                VaccPremDuringCurrentTimestep::Array{Bool,1},
                                VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                VaccEff::Float64,
                                PremVaccStage::Array{Int64,1},
                                prem_time_to_inoculation::Array{Float64,1})

For a single premises, check if vaccination occurs. If yes, cattle are vaccinated.

Inputs:
- `prem_idx::Int64`: ID of premises to be checked.
- `PremNum::Int64`: Number of premises in landscape.
- `DetectionTime::Float64`: Time from infection until premises notifies that infection is present.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.
- `VaccEff::Float64`: Efficacy of vaccine.
- `PremVaccStage::Array{Int64,1}`: Vector signifying the stage premises undergoes vaccination.
- `prem_time_to_inoculation::Array{Float64,1}`: Time for vaccine to become effective. Entry per premises.

Outputs: None \n
Location: ControlFns.jl
"""
function apply_cattle_only_vacc_fn!(prem_idx::Int64,
                                    PremNum::Int64,
                                    DetectionTime::Float64,
                                    PremStatus::Array{Float64,1},
                                    PremVaccStatus::Array{Int64,1},
                                    SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                    PremSuscept::Array{Float64},
                                    PremTransmiss::Array{Float64},
                                    PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                    PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                    VaccPremDuringCurrentTimestep::Array{Bool,1},
                                    VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                    VaccEff::Float64,
                                    PremVaccStage::Array{Int64,1},
                                    prem_time_to_inoculation::Array{Float64,1})

    if ( (PremVaccStatus[prem_idx]) == 0 &&
         (0<=PremStatus[prem_idx]<DetectionTime) &&
         (PremVaccStage[prem_idx] == 2) )

        # Vaccinate premises, update PremVaccStatus
        PremVaccStatus[prem_idx] = 1

        # Only cattle are vaccinated (first column of livestock type)
        # Update SpeciesGroupVaccStatusByPrem
        SpeciesGroupVaccStatusByPrem[prem_idx,1] = 1

        #Update premises vaccinated during current timestep vector
        VaccPremDuringCurrentTimestep[prem_idx] = 1

        # Vaccine only successful in modifying susceptiblity & transmissiblity immediately IF
        # no delay in vaccine becoming effective post being administered AND
        # premises is susceptible
        if (prem_time_to_inoculation[prem_idx] == 0)

            # Check premises is still susceptible.
            # If susceptible, apply vaccination effect
            if (PremStatus[prem_idx] == 0)
                # Modify premises-level values.
                # Update vaccinated premises susceptibility and transmissibility for cattle.
                PremSuscept_SpeciesBreakdown[prem_idx,1] = PremSuscept_SpeciesBreakdown[prem_idx,1].*(1-VaccEff)
                PremTransmiss_SpeciesBreakdown[prem_idx,1] = PremTransmiss_SpeciesBreakdown[prem_idx,1].*(1-VaccEff)
                PremSuscept[prem_idx] = sum(PremSuscept_SpeciesBreakdown[prem_idx,:])
                PremTransmiss[prem_idx] = sum(PremTransmiss_SpeciesBreakdown[prem_idx,:])

                # Update vaccination becomes effective during current timestep vector
                VaccBecomesEffectiveDuringCurrentTimestep[prem_idx] = 1
            end
        end
    end
    return nothing
end

"""
    deploy_cattle_only_vacc_fn!(prem_idx::Int64,
                                DetectionTime::Float64,
                                PremStatus::Array{Float64,1},
                                PremVaccStatus::Array{Int64,1},
                                SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                VaccPremDuringCurrentTimestep::Array{Bool,1},
                                PremVaccStage::Array{Int64,1})

For a single premises, apply cattle vaccination strategy with time lag until effect.

Inputs:
- `prem_idx::Int64`: ID of premises to be checked.
- `DetectionTime::Float64`: Time from infection until premises notifies that infection is present.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `PremVaccStage::Array{Int64,1}`: Vector signifying the stage premises undergoes vaccination.

Outputs: None \n
Location: ControlFns.jl
"""
# - deploy_cattle_only_vacc_fn! (For a single premises, apply cattle vaccination strategy with time lag until effect)
function deploy_cattle_only_vacc_fn!(prem_idx::Int64,
                                    DetectionTime::Float64,
                                    PremStatus::Array{Float64,1},
                                    PremVaccStatus::Array{Int64,1},
                                    SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                    VaccPremDuringCurrentTimestep::Array{Bool,1},
                                    PremVaccStage::Array{Int64,1})

    if ( (PremVaccStatus[prem_idx]) == 0 &&
         (0<=PremStatus[prem_idx]<DetectionTime) &&
         (PremVaccStage[prem_idx] == 2) )

        # Vaccinate premises, update PremVaccStatus
        PremVaccStatus[prem_idx] = 1

        # Only cattle are vaccinated (first column of livestock type)
        # Update SpeciesGroupVaccStatusByPrem
        SpeciesGroupVaccStatusByPrem[prem_idx,1] = 1

        #Update premises vaccinated during current timestep vector
        VaccPremDuringCurrentTimestep[prem_idx] = 1
    end
    return nothing
end

"""
    activate_cattle_only_vacc_fn!(prem_idx::Int64,
                                    PremStatus::Array{Float64,1},
                                    PremSuscept::Array{Float64},
                                    PremTransmiss::Array{Float64},
                                    PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                    PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                    VaccEff::Float64,
                                    VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1})

For a single premises, apply effect of cattle only vaccination strategy.

Inputs:
- `prem_idx::Int64`: ID of premises to be checked.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `VaccEff::Float64`: Efficacy of vaccine.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.

Outputs: None \n
Location: ControlFns.jl
"""
function activate_cattle_only_vacc_fn!(prem_idx::Int64,
                                    PremStatus::Array{Float64,1},
                                    PremSuscept::Array{Float64},
                                    PremTransmiss::Array{Float64},
                                    PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                    PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                    VaccEff::Float64,
                                    VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1})

    # Vaccine only successful in modifying susceptiblity & transmissiblity if susceptible
    if PremStatus[prem_idx] == 0

        # Update vaccinated premises susceptibility and transmissibility for cattle.
        # Modify premises-level values.
        PremSuscept_SpeciesBreakdown[prem_idx,1] = PremSuscept_SpeciesBreakdown[prem_idx,1].*(1-VaccEff)
        PremTransmiss_SpeciesBreakdown[prem_idx,1] = PremTransmiss_SpeciesBreakdown[prem_idx,1].*(1-VaccEff)
        PremSuscept[prem_idx] = sum(PremSuscept_SpeciesBreakdown[prem_idx,:])
        PremTransmiss[prem_idx] = sum(PremTransmiss_SpeciesBreakdown[prem_idx,:])

        # Update vaccination due to become effective during current timestep vector
        VaccBecomesEffectiveDuringCurrentTimestep[prem_idx] = 1
    end
    return nothing
end

"""
    cull_prem_and_update_time_to_all_livestock_vacc_effective_fn!(PremNum::Int64,
                                                                    RemovalTime::Float64,
                                                                    PremStatus::Array{Float64,1},
                                                                    PremVaccStatus::Array{Int64,1},
                                                                    CullPremDuringCurrentTimestep::Array{Bool,1},
                                                                    VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                                                    prem_time_to_inoculation::Array{Float64,1},
                                                                    PremSuscept::Array{Float64},
                                                                    PremTransmiss::Array{Float64},
                                                                    PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                                                    PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                                                    VaccEff::Float64,
                                                                    delta_t::Float64)

Find and update: (i) premises to be culled; (ii) time to vaccine becoming effective (vaccine administered to all livestock types).

Inputs:
- `PremNum::Int64`: Number of premises in landscape.
- `RemovalTime::Float64`: Time from infection until livestock on premises are culled.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.
- `prem_time_to_inoculation::Array{Float64,1}`: Time for vaccine to become effective. Entry per premises.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `VaccEff::Float64`: Efficacy of vaccine.
- `delta_t::Float64`: Timestep.

Outputs: None \n
Location: ControlFns.jl
"""
function cull_prem_and_update_time_to_all_livestock_vacc_effective_fn!(PremNum::Int64,
                                                                        RemovalTime::Float64,
                                                                        PremStatus::Array{Float64,1},
                                                                        PremVaccStatus::Array{Int64,1},
                                                                        CullPremDuringCurrentTimestep::Array{Bool,1},
                                                                        VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                                                        prem_time_to_inoculation::Array{Float64,1},
                                                                        PremSuscept::Array{Float64},
                                                                        PremTransmiss::Array{Float64},
                                                                        PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                                                        PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                                                        VaccEff::Float64,
                                                                        delta_t::Float64
                                                                        )
    for prem_idx = 1:PremNum
        # If applicable, find and update premises to be culled
        if PremStatus[prem_idx]>RemovalTime
            #Cull IP
            PremStatus[prem_idx] = -1.

            #Update premises culled during current timestep vector
            CullPremDuringCurrentTimestep[prem_idx] = 1
        end

        # If applicable, find and update time to vaccine becoming effective
        if PremVaccStatus[prem_idx] == 1 # Find premises that are vaccinated

            # Check if time to inoculation is above zero.
            if prem_time_to_inoculation[prem_idx] > 0
                # If satisfied, decrease by timestep
                prem_time_to_inoculation[prem_idx] -= delta_t

                # Error check
                # prem_time_to_inoculation[prem_idx] should not decrease below zero
                if prem_time_to_inoculation[prem_idx] < 0
                    error("prem_time_to_inoculation[$prem_idx]: $(prem_time_to_inoculation[prem_idx]). Invalid.")
                end

                # If time to inoculation reaches zero, apply effect of vaccination (if valid)
                if prem_time_to_inoculation[prem_idx] == 0

                    # Modify premises-level values for susceptibility and transmissibility for each livestock type.
                    # Vaccine only successful in modifying susceptiblity & transmissiblity if susceptible
                    activate_all_livestock_vacc_fn!(prem_idx,
                                                    PremStatus,
                                                    PremSuscept,
                                                    PremTransmiss,
                                                    PremSuscept_SpeciesBreakdown,
                                                    PremTransmiss_SpeciesBreakdown,
                                                    VaccEff,
                                                    VaccBecomesEffectiveDuringCurrentTimestep)
                end
            end
        end
    end
end


"""
    apply_all_livestock_vacc_fn!(prem_idx::Int64,
                                    PremNum::Int64,
                                    DetectionTime::Float64,
                                    PremStatus::Array{Float64,1},
                                    PremVaccStatus::Array{Int64,1},
                                    SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                    PremSuscept::Array{Float64},
                                    PremTransmiss::Array{Float64},
                                    PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                    PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                    VaccPremDuringCurrentTimestep::Array{Bool,1},
                                    VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                    VaccEff::Float64,
                                    PremVaccStage::Array{Int64,1},
                                    prem_time_to_inoculation::Array{Float64,1})

For a single premises, check if vaccination occurs. If yes, all livestock are vaccinated.

Inputs:
- `prem_idx::Int64`: ID of premises to be checked.
- `PremNum::Int64`: Number of premises in landscape.
- `DetectionTime::Float64`: Time from infection until premises notifies that infection is present.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.
- `VaccEff::Float64`: Efficacy of vaccine.
- `PremVaccStage::Array{Int64,1}`: Vector signifying the stage premises undergoes vaccination.
- `prem_time_to_inoculation::Array{Float64,1}`: Time for vaccine to become effective. Entry per premises.

Outputs: None \n
Location: ControlFns.jl
"""
function apply_all_livestock_vacc_fn!(prem_idx::Int64,
                                        PremNum::Int64,
                                        DetectionTime::Float64,
                                        PremStatus::Array{Float64,1},
                                        PremVaccStatus::Array{Int64,1},
                                        SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                        PremSuscept::Array{Float64},
                                        PremTransmiss::Array{Float64},
                                        PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                        PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                        VaccPremDuringCurrentTimestep::Array{Bool,1},
                                        VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1},
                                        VaccEff::Float64,
                                        PremVaccStage::Array{Int64,1},
                                        prem_time_to_inoculation::Array{Float64,1})

    if ( (PremVaccStatus[prem_idx]) == 0 &&
         (0<=PremStatus[prem_idx]<DetectionTime) &&
         (PremVaccStage[prem_idx] == 2) )

        # Vaccinate premises, update PremVaccStatus
        PremVaccStatus[prem_idx] = 1

        # All livestock types vaccinated
        # Update SpeciesGroupVaccStatusByPrem
        SpeciesGroupVaccStatusByPrem[prem_idx,:] .= 1

        #Update premises vaccinated during current timestep vector
        VaccPremDuringCurrentTimestep[prem_idx] = 1

        # Vaccine only successful in modifying susceptiblity & transmissiblity immediately IF
        # no delay in vaccine becoming effective post being administered AND
        # premises is susceptible
        if (prem_time_to_inoculation[prem_idx] == 0)

            # Check premises is still susceptible.
            # If susceptible, apply vaccination effect
            if (PremStatus[prem_idx] == 0)
                # Modify premises-level values.
                # Update vaccinated premises susceptibility and transmissibility for all livestock types
                PremSuscept_SpeciesBreakdown[prem_idx,:] = PremSuscept_SpeciesBreakdown[prem_idx,:].*(1-VaccEff)
                PremTransmiss_SpeciesBreakdown[prem_idx,:] = PremTransmiss_SpeciesBreakdown[prem_idx,:].*(1-VaccEff)
                PremSuscept[prem_idx] = PremSuscept[prem_idx].*(1-VaccEff)
                PremTransmiss[prem_idx] = PremTransmiss[prem_idx].*(1-VaccEff)

                # Update vaccination becomes effective during current timestep vector
                VaccBecomesEffectiveDuringCurrentTimestep[prem_idx] = 1
            end
        end
    end
    return nothing
end


"""
    deploy_all_livestock_vacc_fn!(prem_idx::Int64,
                                    DetectionTime::Float64,
                                    PremStatus::Array{Float64,1},
                                    PremVaccStatus::Array{Int64,1},
                                    SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                    VaccPremDuringCurrentTimestep::Array{Bool,1},
                                    PremVaccStage::Array{Int64,1})

For a single premises, apply livestock vaccination strategy with time lag until effect.

Inputs:
- `prem_idx::Int64`: ID of premises to be checked.
- `DetectionTime::Float64`: Time from infection until premises notifies that infection is present.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `PremVaccStage::Array{Int64,1}`: Vector signifying the stage premises undergoes vaccination.

Outputs: None \n
Location: ControlFns.jl
"""
function deploy_all_livestock_vacc_fn!(prem_idx::Int64,
                                        DetectionTime::Float64,
                                        PremStatus::Array{Float64,1},
                                        PremVaccStatus::Array{Int64,1},
                                        SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                        VaccPremDuringCurrentTimestep::Array{Bool,1},
                                        PremVaccStage::Array{Int64,1})


    if ( (PremVaccStatus[prem_idx]) == 0 &&
         (0<=PremStatus[prem_idx]<DetectionTime) &&
         (PremVaccStage[prem_idx] == 2) )

        # Vaccinate premises, update PremVaccStatus
        PremVaccStatus[prem_idx] = 1

        # All livestock types vaccinated
        # Update SpeciesGroupVaccStatusByPrem
        SpeciesGroupVaccStatusByPrem[prem_idx,:] .= 1

        #Update premises vaccinated during current timestep vector
        VaccPremDuringCurrentTimestep[prem_idx] = 1
    end
    return nothing
end

"""
    activate_all_livestock_vacc_fn!(prem_idx::Int64,
                                    PremStatus::Array{Float64,1},
                                    PremSuscept::Array{Float64},
                                    PremTransmiss::Array{Float64},
                                    PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                    PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                    VaccEff::Float64,
                                    VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1})

For a single premises, apply effect of livestock vaccination strategy.

Inputs:
- `prem_idx::Int64`: ID of premises to be checked.
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `VaccEff::Float64`: Efficacy of vaccine.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.

Outputs: None \n
Location: ControlFns.jl
"""
function activate_all_livestock_vacc_fn!(prem_idx::Int64,
                                            PremStatus::Array{Float64,1},
                                            PremSuscept::Array{Float64},
                                            PremTransmiss::Array{Float64},
                                            PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                            PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                            VaccEff::Float64,
                                            VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1})

    # Vaccine only successful in modifying susceptiblity & transmissiblity if susceptible
    if PremStatus[prem_idx] == 0

        # Update vaccinated premises susceptibility and transmissibility for cattle.
        # Modify premises-level values.
        PremSuscept_SpeciesBreakdown[prem_idx,:] = PremSuscept_SpeciesBreakdown[prem_idx,:].*(1-VaccEff)
        PremTransmiss_SpeciesBreakdown[prem_idx,:] = PremTransmiss_SpeciesBreakdown[prem_idx,:].*(1-VaccEff)
        # PremSuscept[prem_idx] = sum(PremSuscept_SpeciesBreakdown[prem_idx,:])
        # PremTransmiss[prem_idx] = sum(PremTransmiss_SpeciesBreakdown[prem_idx,:])
        PremSuscept[prem_idx] = PremSuscept[prem_idx].*(1-VaccEff)
        PremTransmiss[prem_idx] = PremTransmiss[prem_idx].*(1-VaccEff)

        # Update vaccination due to become effective during current timestep vector
        VaccBecomesEffectiveDuringCurrentTimestep[prem_idx] = 1
    end
    return nothing
end
