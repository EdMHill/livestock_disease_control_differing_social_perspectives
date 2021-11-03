#=
Purpose:
This file houses functions to perform single replicate of spatial simulation using the Sellke construction

Date: 3rd November 2021
=#

#-------------------------------------------------------------------------------
# UPDATE PREMISES STATUS FOLLOWING INFECTION UPDATES & BEFORE CONTROL CARRIED OUT
#-------------------------------------------------------------------------------
"""
    UpdatePremStatus!(PremStatus::Array{Float64,1},delta_t::Float64,PremHasHadInfFlag::Array{Int64,1})

Construct a collection of Binomial RNGs (random number generators).

Inputs:
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `delta_t::Float64`: Timestep increment.
- `PremHasHadInfFlag::Array{Int64,1}`: Premises-level record of a premises ever having had infection (previously or currently). Entry per premises. Never infected: 0. Has been infected: 1.

Outputs: None \n
Location: RunSellkeOutbreakReplicateFns.jl
"""
function UpdatePremStatus!(PremStatus::Array{Float64,1},
                            delta_t::Float64,
                            PremHasHadInfFlag::Array{Int64,1})

    for PremIdxItr = 1:length(PremStatus)
        if PremStatus[PremIdxItr]>0 #Premises previously infected. Update status
            PremStatus[PremIdxItr] = PremStatus[PremIdxItr] + delta_t
        elseif PremStatus[PremIdxItr]==0 #Premises susceptible. Add infection flag value. Those infected on current timestep will update
            PremStatus[PremIdxItr] = PremStatus[PremIdxItr] + (PremHasHadInfFlag[PremIdxItr]*delta_t)
        end
    end
    return nothing
end

#-------------------------------------------------------------------------------
# PERFORM SINGLE TIMESTEP OF OUTBREAK
#-------------------------------------------------------------------------------
"""
    IterateOutbreak_TestSellkeFMDsimn!(EpiParamVals::Array{Float64,1},
                                        ControlParamVals::Array{Any,1},
                                        LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                        PremNum::Int64,
                                        PremLoc_AllVals::Array{Float64,2},
                                        PremLoc_xVals::Array{Float64},
                                        PremLoc_yVals::Array{Float64},
                                        PremSuscept::Array{Float64},
                                        PremTransmiss::Array{Float64},
                                        PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                        PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                        PremStatus::Array{Float64,1},
                                        PremHasHadInfFlag::Array{Int64,1},
                                        notified_prem_timeseries::Array{Int64,1},
                                        PremVaccStatus::Array{Int64,1},
                                        SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                        PremInfectiousFlagVec::Array{Int64,1},
                                        CoordType::Int64,
                                        KernelLookUpVec::Array{Float64,1},
                                        RunControlsFn::Function,
                                        delta_t::Float64,
                                        prem_suscep_remaining_sellke::Array{Float64,1},
                                        precalc_dist_and_kernel_array_flag::Bool,
                                        kernel_val_array::Array{Float64,2})

Perform a single timestep of the spatial outbreak simulation using the Sellke construction.

Inputs:
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)]
- `ControlParamVals::Array{Any,1}`: Variables related to implementing control measures
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control
- `PremNum::Int64`: Number of premises in landscape.
- `PremLoc_AllVals::Array{Float64,2}`: Coordinates for each premises.
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremHasHadInfFlag::Array{Int64,1}`: Premises-level record of a premises ever having had infection (previously or currently). Entry per premises. Never infected: 0. Has been infected: 1.
- `notified_prem_timeseries::Array{Int64,1}`: Per timestep, the number of premises that reported infection
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `PremInfectiousFlagVec::Array{Int64,1}`: Indicator of whether node is currently infected or not. As vector it has an entry per premises. Can have as array if infection over premises is not complete. Row per premises, column per group.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `KernelLookUpVec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `RunControlsFn::Function`: Function containing control measures to be enacted.
- `delta_t::Float64`: Timestep between each iteration.
- `prem_suscep_remaining_sellke::Array{Float64,1}`: Susceptibility/Resistance remaining for each premises.
- `precalc_dist_and_kernel_array_flag::Bool`: Indicator that if true, construct array of premises-to-premises kernel quantities before running main outbreak function.
- `kernel_val_array::Array{Float64,2}`:  Distance kernel evaluated for each premises-to-premises combination.

Outputs:
- `CullPremDuringCurrentTimestep::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on premises during current timestep, 0 otherwise.
- `VaccPremDuringCurrentTimestep::Array{Bool,1}`: Indicator vector. 1 if premises culled during current timestep. 0 otherwise.
- `VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.

Location: RunSellkeOutbreakReplicateFns.jl
"""
function IterateOutbreak_TestSellkeFMDsimn!(EpiParamVals::Array{Float64,1},
                                            ControlParamVals::Array{Any,1},
                                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                            PremNum::Int64,
                                            PremLoc_AllVals::Array{Float64,2},
                                            PremLoc_xVals::Array{Float64},
                                            PremLoc_yVals::Array{Float64},
                                            PremSuscept::Array{Float64},
                                            PremTransmiss::Array{Float64},
                                            PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                            PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                            PremStatus::Array{Float64,1},
                                            PremHasHadInfFlag::Array{Int64,1},
                                            notified_prem_timeseries::Array{Int64,1},
                                            PremVaccStatus::Array{Int64,1},
                                            SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                            PremInfectiousFlagVec::Array{Int64,1},
                                            CoordType::Int64,
                                            KernelLookUpVec::Array{Float64,1},
                                            RunControlsFn::Function,
                                            delta_t::Float64,
                                            prem_suscep_remaining_sellke::Array{Float64,1},
                                            precalc_dist_and_kernel_array_flag::Bool,
                                            kernel_val_array::Array{Float64,2})

    #---------------------------------------------------------------------------
    # RUN TRANSMISSION PROCESS
    #---------------------------------------------------------------------------

    #Find infectious premises. Remove any repeat entries
    NumInfPrem = sum(PremInfectiousFlagVec.==1)

    # Check if any infections
    if NumInfPrem == 0 #If no infectious premises, update status for those exposed
        #Update status for those previously infected
        for PremIdxItr = 1:PremNum
            if PremStatus[PremIdxItr]>0 #Premises previously infected. Update status
                PremStatus[PremIdxItr] = PremStatus[PremIdxItr] + delta_t
            end
        end
    else
        # If infectious premises present,
        #   - get list of susceptible premises
        #   - iterate over each susceptible premises, calculate FOI from each infectious premises (and kernal value if needed)
        #   - modify remaining suscepbility/resistance for susceptible premises
        #   - After all FOI applied, check for resistances in susceptible premises that have gone below 0. Update status
        #   - End of timestep, pre-control: Update status for premises
        #   - Initiate controls at end of timestep

        # Get list of susceptible premises
        suscep_prem_ID_vec = findall(PremStatus .== 0) #Get IDs of premises that are susceptible
        n_suscep_prem = length(suscep_prem_ID_vec) #Number of susceptible premises

        # Get list of infectious premises
        infectious_prem_ID_vec = findall(PremInfectiousFlagVec.==1)
        n_infectious_prem = length(infectious_prem_ID_vec)

        # Iterate over each susceptible premises, calculate FOI from each infectious premises (and kernal value if needed)
        for suscep_prem_itr = 1:n_suscep_prem

            # Get ID of susceptible premises
            suscep_prem_ID = suscep_prem_ID_vec[suscep_prem_itr]

            # Calculate kernel value if NOT precalculated
            if precalc_dist_and_kernel_array_flag == false
                # Assign location info to variable if premises-to-premises distances
                # are being calculated on the fly
                suscep_prem_xLoc = PremLoc_xVals[suscep_prem_ID]
                suscep_prem_yLoc = PremLoc_yVals[suscep_prem_ID]
            end

            # Get FOI against the susceptible premises provided by each infectious premises
            for infec_prem_itr = 1:n_infectious_prem

                # Get ID of infectious premises
                infectious_prem_ID = infectious_prem_ID_vec[infec_prem_itr]

                # Calculate kernel value if NOT precalculated
                if precalc_dist_and_kernel_array_flag == false

                    # Assign location info to variable if premises-to-premises distances
                    # are being calculated on the fly
                    infectious_prem_xLoc = PremLoc_xVals[infectious_prem_ID]
                    infectious_prem_yLoc = PremLoc_yVals[infectious_prem_ID]

                    #Calculate distance between the two points
                    if CoordType == 1 #Cartesian co-ords (metres)
                        d =  eucl_distance(infectious_prem_xLoc,
                                            infectious_prem_yLoc,
                                            suscep_prem_xLoc,
                                            suscep_prem_yLoc)
                    elseif CoordType == 2 #Cartesian co-ords (metres)
                        d = eucl_distance_ConvertToMetres(infectious_prem_xLoc,
                                                            infectious_prem_yLoc,
                                                            suscep_prem_xLoc,
                                                            suscep_prem_yLoc)
                    elseif CoordType == 3 #Lat/Long co-ords
                        d = GreatCircleDistance(infectious_prem_yLoc, infectious_prem_xLoc,  #lat1, lon1
                                                suscep_prem_yLoc, suscep_prem_xLoc) #lat2, lon2
                    end

                    # Get distance between the two premises
                    distIdx = ReturnDistIdxForKernel(d)

                    # Apply the kernel to that distance
                    kernel_value = KernelLookUpVec[distIdx]
                else
                    kernel_value = kernel_val_array[suscep_prem_ID,infectious_prem_ID]
                end

                # Calculate force of infection
                FOIrate = PremTransmiss[infectious_prem_ID]*PremSuscept[suscep_prem_ID]*kernel_value

                # Modify remaining suscepbility/resistance for susceptible premises
                prem_suscep_remaining_sellke[suscep_prem_ID] -= oneMinusExp(-FOIrate*delta_t)
            end

            # After all FOI applied, check for resistances in susceptible premises that have gone below 0. Update status
            if prem_suscep_remaining_sellke[suscep_prem_ID] < 0
                PremHasHadInfFlag[suscep_prem_ID] = 1
            end
        end

        #End of timestep, pre-control: Update status for premises
        UpdatePremStatus!(PremStatus,delta_t,PremHasHadInfFlag)
    end

    #---------------------------------------------------------------------------
    # RUN CONTROLS
    #---------------------------------------------------------------------------

    #Initiate controls at end of timestep
    CullPremDuringCurrentTimestep = zeros(Bool,PremNum)
    if (RunControlsFn == CullIPsAndRingVaccFn!)  #Vaccine based strategy, all livestock types targeted

        #Initialise flag vector for vaccination being administered
        # during current timestep
        VaccPremDuringCurrentTimestep = zeros(Bool,PremNum)

        #Initialise flag vector for vaccination becoming effective current timestep
        VaccBecomesEffectiveDuringCurrentTimestep = zeros(Bool,PremNum)

        #Run control function
        RunControlsFn(PremNum,
                        PremStatus,
                        PremVaccStatus,
                        SpeciesGroupVaccStatusByPrem,
                        PremSuscept,
                        PremTransmiss,
                        PremLoc_xVals,
                        PremLoc_yVals,
                        EpiParamVals,
                        ControlParamVals,
                        LinkedPremToControlIdxs,
                        CoordType,
                        CullPremDuringCurrentTimestep,
                        VaccPremDuringCurrentTimestep,
                        VaccBecomesEffectiveDuringCurrentTimestep)
        return CullPremDuringCurrentTimestep::Array{Bool,1},
                VaccPremDuringCurrentTimestep::Array{Bool,1},
                VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}
    elseif ( (RunControlsFn == CullIPsAndRingVaccCattleOnlyFn!) ||
            (RunControlsFn == CullIPsAndDistanceBasedVaccFn!) ||
            (RunControlsFn == CullIPsAndDistanceBasedCattleVaccFn!) )
        # Vaccine based strategy, specific livestock types targeted

        # Initialise flag vector for vaccination during current timestep
        VaccPremDuringCurrentTimestep = zeros(Bool,PremNum)

        #Initialise flag vector for vaccination becoming effective current timestep
        VaccBecomesEffectiveDuringCurrentTimestep = zeros(Bool,PremNum)

        # Run control function
        RunControlsFn(PremNum,
                        PremStatus,
                        PremVaccStatus,
                        SpeciesGroupVaccStatusByPrem,
                        PremSuscept,
                        PremTransmiss,
                        PremSuscept_SpeciesBreakdown,
                        PremTransmiss_SpeciesBreakdown,
                        PremLoc_xVals,
                        PremLoc_yVals,
                        EpiParamVals,
                        ControlParamVals,
                        LinkedPremToControlIdxs,
                        CoordType,
                        CullPremDuringCurrentTimestep,
                        VaccPremDuringCurrentTimestep,
                        VaccBecomesEffectiveDuringCurrentTimestep,
                        delta_t)
                        return CullPremDuringCurrentTimestep::Array{Bool,1},
                                VaccPremDuringCurrentTimestep::Array{Bool,1},
                                VaccBecomesEffectiveDuringCurrentTimestep::Array{Bool,1}
    else
        RunControlsFn(PremNum,
                        PremStatus,
                        PremLoc_xVals,
                        PremLoc_yVals,
                        EpiParamVals,
                        ControlParamVals,
                        LinkedPremToControlIdxs,
                        CoordType,
                        CullPremDuringCurrentTimestep)
        return CullPremDuringCurrentTimestep::Array{Bool,1}
    end
end

#-------------------------------------------------------------------------------
# PERFORM COMPLETE OUTBREAK REPLICATE
#-------------------------------------------------------------------------------
"""
    RunTestSellkeFMDsimn(BatchID::String,
                            ItrIdx::Int64,
                            EpiParamVals::Array{Float64,1},
                            ControlParamVals::Array{Any,1},
                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                            CoordType::Int64,
                            PremStatus::Array{Float64,1},
                            PremHasHadInfFlag::Array{Int64,1},
                            PremVaccStatus::Array{Int64,1},
                            SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                            EventIndicatorArray::Array{Int64,2},
                            EventTimeArray::Array{Float64,2},
                            PremLoc_AllVals::Array{Float64,2},
                            PremLoc_xVals::Array{Float64,1},
                            PremLoc_yVals::Array{Float64,1},
                            PremLivestockData::Array{Int64},
                            PremSuscept::Array{Float64},
                            PremTransmiss::Array{Float64},
                            PremSuscept_SpeciesBreakdown::Array{Float64,2},
                            PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                            OutputLevels, #Either Int64, value 0, or an Array{Int64,1}.
                            MaxTime::Float64,
                            TimeStepVal::Float64,
                            KernelLookUpVec::Array{Float64,1},
                            RunControlsFn::Function,
                            ReplicateFilePrefix::String,
                            OutputFileObjs::Array{IOStream,1},
                            IterateOutbreakFn::Function,
                            prem_suscep_remaining_sellke::Array{Float64,1},
                            precalc_dist_and_kernel_array_flag::Bool,
                            kernel_val_array::Array{Float64,2})

Perform a complete replicate of the spatial outbreak simulation using the Sellke construction.

Inputs:
- `BatchID::String`: Used as prefix for output files.
- `ItrIdx::Int64`: ID for replicate
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)]
- `ControlParamVals::Array{Any,1}`: Variables related to implementing control measures
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `PremStatus::Array{Float64,1}`: Premises-level disease status. Entry per premises. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `PremHasHadInfFlag::Array{Int64,1}`: Premises-level record of a premises ever having had infection (previously or currently). Entry per premises. Never infected: 0. Has been infected: 1.
- `PremVaccStatus::Array{Int64,1}`: Premises-level vaccination status. Entry per premises. Unvaccinated: 0. Vaccinated: 1.
- `SpeciesGroupVaccStatusByPrem::Array{Int64,2}`: At each premises, livestock group vaccination status. Unvaccinated: 0. Vaccinated: 1.
- `EventIndicatorArray::Array{Int64,2}`: Flag recording if event occurs. Row per premises. First column gives premises ID.
- `EventTimeArray::Array{Float64,2}`: Array recording time the event occurred at. Row per premises. First column gives premises ID.
- `PremLoc_AllVals::Array{Float64,2}`: Coordinates for each premises.
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `PremLivestockData::Array{Int64}`: Number of each livestock type per premises. Row per premises, column per animal.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`: Per premises, the susceptibility values for each speices/livestock type (row per premises, column per species).
- `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, the transmissibility values for each speices/livestock type (row per premises, column per species).
- `OutputLevels`: Will output replicate results when the cumulative number of infected reaches these levels. Will also exit after the last output. Set to 0 if want to output every timestep.
- `MaxTime::Float64`: Timeframe simulation is run over
- `TimeStepVal::Float64`: Timestep per iteration
- `KernelLookUpVec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `RunControlsFn::Function`: Enact controls as specified within the given function.
- `ReplicateFilePrefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `OutputFileObjs::Array{IOStream,1}`: Filename identifiers written to by all replicates.
- `IterateOutbreakFn::Function`: Function to perform disease transitions and control implementation per timestep
- `prem_suscep_remaining_sellke::Array{Float64,1}`: Susceptibility/Resistance remaining for each premises.
- `precalc_dist_and_kernel_array_flag::Bool`: Indicator that if true, construct array of premises-to-premises kernel quantities before running main outbreak function.
- `kernel_val_array::Array{Float64,2}`:  Distance kernel evaluated for each premises-to-premises combination.

Outputs: None \n
Location: RunSellkeOutbreakReplicateFns.jl
"""
function RunTestSellkeFMDsimn(BatchID::String,
                                ItrIdx::Int64,
                                EpiParamVals::Array{Float64,1},
                                ControlParamVals::Array{Any,1},
                                LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                CoordType::Int64,
                                PremStatus::Array{Float64,1},
                                PremHasHadInfFlag::Array{Int64,1},
                                PremVaccStatus::Array{Int64,1},
                                SpeciesGroupVaccStatusByPrem::Array{Int64,2},
                                EventIndicatorArray::Array{Int64,2},
                                EventTimeArray::Array{Float64,2},
                                PremLoc_AllVals::Array{Float64,2},
                                PremLoc_xVals::Array{Float64,1},
                                PremLoc_yVals::Array{Float64,1},
                                PremLivestockData::Array{Int64},
                                PremSuscept::Array{Float64},
                                PremTransmiss::Array{Float64},
                                PremSuscept_SpeciesBreakdown::Array{Float64,2},
                                PremTransmiss_SpeciesBreakdown::Array{Float64,2},
                                OutputLevels, #Either Int64, value 0, or an Array{Int64,1}.
                                MaxTime::Float64,
                                TimeStepVal::Float64,
                                KernelLookUpVec::Array{Float64,1},
                                RunControlsFn::Function,
                                ReplicateFilePrefix::String,
                                OutputFileObjs::Array{IOStream,1},
                                IterateOutbreakFn::Function,
                                prem_suscep_remaining_sellke::Array{Float64,1},
                                precalc_dist_and_kernel_array_flag::Bool,
                                kernel_val_array::Array{Float64,2})

    #---------------------------------------------------------------------------
    #Unpack Epidemioloigcal parameters
    #---------------------------------------------------------------------------
    IncubationTime = EpiParamVals[1]::Float64
    DetectionTime = EpiParamVals[2]::Float64
    RemovalTime = EpiParamVals[3]::Float64

    #---------------------------------------------------------------------------
    #Intiailise variables to be used throughout script
    #---------------------------------------------------------------------------

    #Get number of premises in use
    PremNum = size(PremStatus,1)

    #Get number of livestock species in use
    LivestockTypes = size(PremLivestockData,2)

    #---------------------------------------------------------------------------
    #Open files to store outputs, for specific replicate
    #---------------------------------------------------------------------------
    EventIndicatorArrayFile = string(ReplicateFilePrefix,"_EventArrays/EventIndicatorArray_BatchID$(BatchID)_Replicate$(ItrIdx).txt")
    EventTimeArrayFile = string(ReplicateFilePrefix,"_EventArrays/EventTimeArray_BatchID$(BatchID)_Replicate$(ItrIdx).txt")

    PremPerDiseaseStatePerTimestepFile = string(ReplicateFilePrefix,"_PerTimestep/PremPerDiseaseStatePerTimestep_BatchID$(BatchID)_Replicate$(ItrIdx).txt")
    CumulativeCulledPerTimestepFile = string(ReplicateFilePrefix,"_PerTimestep/CumulativeCulledPerTimestep_BatchID$(BatchID)_Replicate$(ItrIdx).txt")
    CumulVaccPerTimestepFile = string(ReplicateFilePrefix,"_PerTimestep/CumulVaccPerTimestepFile_BatchID$(BatchID)_Replicate$(ItrIdx).txt")
    CumulCasesAnimalsPerTimestepFile = string(ReplicateFilePrefix,"_PerTimestep/CumulCasesAnimalsPerTimestepFile_BatchID$(BatchID)_Replicate$(ItrIdx).txt")

    #---------------------------------------------------------------------------
    #Initialise simulation tracking variables
    #---------------------------------------------------------------------------
    t = 0.
    SimnItr = 1

    #---------------------------------------------------------------------------v
    #Initialise compartment status variables
    #---------------------------------------------------------------------------
    ReqItr = convert(Int64,ceil(MaxTime/TimeStepVal) + 1)  #Add 1 as also storing value at intial time, then every TimeStepVal until MaxTime exceeded
    println("ReqItr: $ReqItr")
    S = zeros(Int64, ReqItr); E = zeros(Int64, ReqItr); I = zeros(Int64, ReqItr);
    R = zeros(Int64, ReqItr); R2 = zeros(Int64, ReqItr); Culled = zeros(Int64, ReqItr);

    # Premises level storage vectors
    VaccPrem = zeros(Int64, ReqItr); CumulCasesPrem = zeros(Int64, ReqItr);
    CumulCasesWithoutControlPrem = zeros(Int64, ReqItr);
    CumulCasesWithControlPrem = zeros(Int64, ReqItr);

    # Animal level storage vectors
    CullAnimals = zeros(Int64,ReqItr,LivestockTypes)
    VaccAnimals = zeros(Int64,ReqItr,LivestockTypes)
    CumulCasesAnimals = zeros(Int64,ReqItr,LivestockTypes)
    CumulCasesWithoutControlAnimals = zeros(Int64,ReqItr,LivestockTypes)
    CumulCasesWithControlAnimals = zeros(Int64,ReqItr,LivestockTypes)

    # Assign IDs of initially infected units to variable
    Infected = findall(IncubationTime .<PremStatus .<=RemovalTime)

    #---------------------------------------------------------------------------
    #Assign first entry to compartment status variables
    #---------------------------------------------------------------------------
    #Equivalent to:
        # S[SimnItr] = sum(PremStatus.==0)
        # E[SimnItr] = sum(0 .< PremStatus .<=IncubationTime)
        # I[SimnItr] = sum(IncubationTime .<PremStatus .<=DetectionTime)
        # R[SimnItr] = sum(PremStatus.==(DetectionTime+1))
        # R2[SimnItr] = sum(PremStatus.>DetectionTime)
        # Culled[SimnItr] = sum(PremStatus.<0)
        # CullAnimals[SimnItr,:] = zeros(Int64,size(PremLivestockData,2))

    #Initialise temporary storage vectors
    LogicTempIntVec = Array{Int64,1}(undef,PremNum)
    LogicTempIntVec2 = Array{Int64,1}(undef,PremNum)
    LogicTempIntArray = Array{Int64,2}(undef,PremNum,LivestockTypes)
    LogicTempVec = BitArray{1}(undef,PremNum)

    #Susceptibles
    PremSusFlagVec = PremStatus .== 0 #Initialise susceptibility flag values
    S[SimnItr] = sum(PremSusFlagVec) #Susceptible status

    #Exposed status. Two conditions to check
    LogicTempIntVec .= PremStatus .<=IncubationTime
    LogicTempIntVec2 .= 0 .< PremStatus #Status must be above 0
    LogicTempIntVec .= LogicTempIntVec.*LogicTempIntVec2 #Take product of two conditions. Those that satisfy both will return value 1
    E[SimnItr]  = sum(LogicTempIntVec)

    #Infectious status. Two conditions to check
    LogicTempIntVec .= PremStatus .<=DetectionTime
    LogicTempIntVec2 .= IncubationTime .<PremStatus #Status must be above 0
    LogicTempIntVec .= LogicTempIntVec.*LogicTempIntVec2 #Take product of two conditions. Those that satisfy both will return value 1
    I[SimnItr]  = sum(LogicTempIntVec)

    #Reported during current time step
    LogicTempIntVec .= PremStatus.==(DetectionTime+1)
    R[SimnItr]  = sum(LogicTempIntVec)

    #All reported premises
    LogicTempIntVec .= PremStatus.>DetectionTime
    R2[SimnItr]  = sum(LogicTempIntVec)

    #Culled
    LogicTempIntVec .= PremStatus.<0
    Culled[SimnItr]  = sum(LogicTempIntVec)

    #Culled animal count
    LogicTempVec .= LogicTempIntVec.==1
    CullAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

    #Vaccinated premises
    LogicTempIntVec .= PremVaccStatus.==1
    VaccPrem[SimnItr]  = sum(LogicTempIntVec)

    #Vaccinated animal count
    LogicTempIntArray .= SpeciesGroupVaccStatusByPrem.*PremLivestockData
    VaccAnimals[SimnItr,:] = sum(LogicTempIntArray,dims=1)

    #Cumulative cases. Two conditions to check
    CumulCasesPrem[SimnItr] = sum(PremHasHadInfFlag)
    LogicTempVec .= PremHasHadInfFlag.==1
    CumulCasesAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

    # Cumulative cases, control NOT applied at holding over entire simulation
    LogicTempIntVec .= PremVaccStatus.==0 # Control has not been applied
    LogicTempVec .= PremHasHadInfFlag.==1
    LogicTempIntVec .= LogicTempIntVec.*LogicTempVec #Take product of two conditions. Those that satisfy both will return value 1
    CumulCasesWithoutControlPrem[SimnItr] = sum(LogicTempIntVec)
    LogicTempVec .= LogicTempIntVec.==1
    CumulCasesWithoutControlAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

    # Cumulative cases, control was applied at holding by end of simulation
    LogicTempIntVec .= PremVaccStatus.*PremHasHadInfFlag
    CumulCasesWithControlPrem[SimnItr] = sum(LogicTempIntVec)
    LogicTempVec .= LogicTempIntVec.==1
    CumulCasesWithControlAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

    #---------------------------------------------------------------------------
    #Initialise infectious flag values
    #---------------------------------------------------------------------------
    PremInfectiousFlagVec = zeros(Int64,PremNum) #Entry per premises
    PremInfectiousFlagVec[Infected] .= 1  #For those premises infected, amend flag value to 1

    #-----------------------------------------------------------------------
    #Update EventIndicatorArray & EventTimeArray
    #-----------------------------------------------------------------------

    #Seed sites, latently infected at initial timestep
    LogicTempVec .= PremStatus.==TimeStepVal
    EventIndicatorArray[LogicTempVec,2] .= 1
    EventTimeArray[LogicTempVec,2] .= t

    # For premises vaccinated at start of simulation, update event arrays
    if ( (RunControlsFn == CullIPsAndDistanceBasedVaccFn!) ||
         (RunControlsFn == CullIPsAndDistanceBasedCattleVaccFn!))

        # Record that vaccine has been given
        LogicTempVec .= PremVaccStatus.==1
        EventIndicatorArray[LogicTempVec,6] .= 1
        EventTimeArray[LogicTempVec,6] .= t

        # Assumed delay to vaccine becoming effective would be surpassed
        # Record that vaccine is effective
        LogicTempVec .= PremVaccStatus.==1
        EventIndicatorArray[LogicTempVec,7] .= 1
        EventTimeArray[LogicTempVec,7] .= t
    end

    #---------------------------------------------------------------------------
    #Initialise OutputLevelIdx if required
    #---------------------------------------------------------------------------
    if ndims(OutputLevels) == 1
        OutputLevelIdx = 1
        CurrentReqOutputLevel = OutputLevels[OutputLevelIdx]
    elseif ndims(OutputLevels) > 1
        error("Dimension of OutputLevels is $OutputLevels. Should have dimension 0 or 1.")
    end

    #---------------------------------------------------------------------------
    #Enter main iteration
    #---------------------------------------------------------------------------
    SimnItr = SimnItr + 1
    IterateFlag = 1 #Flag variable that when altered will lead exit from loop
    while (t < MaxTime && IterateFlag == 1)
        #Call event update function (outputs depend on chosen control fn)
        if ( (RunControlsFn == CullIPsAndDistanceBasedVaccFn!) ||
            (RunControlsFn == CullIPsAndDistanceBasedCattleVaccFn!) ||
            (RunControlsFn == CullIPsAndRingVaccFn!) ||
            (RunControlsFn == CullIPsAndRingVaccCattleOnlyFn!) ||
            (RunControlsFn == CullIPsAndCumulPremCaseBasedVaccFn!) ||
            (RunControlsFn == CullIPsAndPremPrevalenceBasedVaccFn!) )
            #Vaccine based strategy

            CullPremDuringCurrentTimestep,
            VaccPremDuringCurrentTimestep,
            VaccBecomesEffectiveDuringCurrentTimestep = IterateOutbreakFn(EpiParamVals,
                                                                            ControlParamVals,
                                                                            LinkedPremToControlIdxs,
                                                                            PremNum,
                                                                            PremLoc_AllVals,
                                                                            PremLoc_xVals,
                                                                            PremLoc_yVals,
                                                                            PremSuscept,
                                                                            PremTransmiss,
                                                                            PremSuscept_SpeciesBreakdown,
                                                                            PremTransmiss_SpeciesBreakdown,
                                                                            PremStatus,
                                                                            PremHasHadInfFlag,
                                                                            R,
                                                                            PremVaccStatus,
                                                                            SpeciesGroupVaccStatusByPrem,
                                                                            PremInfectiousFlagVec,
                                                                            CoordType,
                                                                            KernelLookUpVec,
                                                                            RunControlsFn,
                                                                            TimeStepVal,
                                                                            prem_suscep_remaining_sellke,
                                                                            precalc_dist_and_kernel_array_flag,
                                                                            kernel_val_array)
        else # Strategy not involving vaccination
            CullPremDuringCurrentTimestep = IterateOutbreakFn(EpiParamVals,
                                                                ControlParamVals,
                                                                LinkedPremToControlIdxs,
                                                                PremNum,
                                                                PremLoc_AllVals,
                                                                PremLoc_xVals,
                                                                PremLoc_yVals,
                                                                PremSuscept,
                                                                PremTransmiss,
                                                                PremSuscept_SpeciesBreakdown,
                                                                PremTransmiss_SpeciesBreakdown,
                                                                PremStatus,
                                                                PremHasHadInfFlag,
                                                                R,
                                                                PremVaccStatus,
                                                                SpeciesGroupVaccStatusByPrem,
                                                                PremInfectiousFlagVec,
                                                                CoordType,
                                                                KernelLookUpVec,
                                                                RunControlsFn,
                                                                TimeStepVal,
                                                                prem_suscep_remaining_sellke,
                                                                precalc_dist_and_kernel_array_flag,
                                                                kernel_val_array)
        end

        #-----------------------------------------------------------------------
        ### Amend compartment based on status
        #-----------------------------------------------------------------------
        # Equivalent to:
            # S[SimnItr] = sum(PremStatus.==0)
            # E[SimnItr] = sum(0 .< PremStatus .<=IncubationTime)
            # I[SimnItr] = sum(IncubationTime .<PremStatus .<=DetectionTime)
            # R[SimnItr] = sum(PremStatus.==(DetectionTime+1))
            # R2[SimnItr] = sum(PremStatus.>DetectionTime)
            # Culled[SimnItr] = sum(PremStatus.<0)

        #Update node type variables, suscept flag values

        #Susceptible status
        PremSusFlagVec .= PremStatus.==0.
        S[SimnItr] = sum(PremSusFlagVec)

        #Exposed status. Two conditions to check
        LogicTempIntVec .= PremStatus .<=IncubationTime
        LogicTempIntVec2 .= 0 .< PremStatus #Status must be above 0
        LogicTempIntVec .= LogicTempIntVec.*LogicTempIntVec2 #Take product of two conditions. Those that satisfy both will return value 1
        E[SimnItr]  = sum(LogicTempIntVec)

        #Infectious status. Two conditions to check
        LogicTempIntVec .= PremStatus .<=DetectionTime
        LogicTempIntVec2 .= IncubationTime .<PremStatus #Status must be above 0
        LogicTempIntVec .= LogicTempIntVec.*LogicTempIntVec2 #Take product of two conditions. Those that satisfy both will return value 1
        I[SimnItr]  = sum(LogicTempIntVec)

        #Reported during current timestep
        LogicTempIntVec .= PremStatus.==(DetectionTime+TimeStepVal)
        R[SimnItr]  = sum(LogicTempIntVec)

        #All reported
        LogicTempIntVec .= PremStatus.>DetectionTime
        R2[SimnItr]  = sum(LogicTempIntVec)

        #Culled
        LogicTempIntVec .= PremStatus.<0
        Culled[SimnItr]  = sum(LogicTempIntVec)

        LogicTempVec .= LogicTempIntVec.==1
        CullAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

        #Vaccinated
        LogicTempIntVec .= PremVaccStatus.==1
        VaccPrem[SimnItr]  = sum(LogicTempIntVec)

        # LogicTempVec .= LogicTempIntVec.==1
        # VaccAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

        #Vaccinated animal count
        LogicTempIntArray .= SpeciesGroupVaccStatusByPrem.*PremLivestockData
        VaccAnimals[SimnItr,:] = sum(LogicTempIntArray,dims=1)

        #Cumulative cases. Two conditions to check
        CumulCasesPrem[SimnItr] = sum(PremHasHadInfFlag)
        LogicTempVec .= PremHasHadInfFlag.==1
        CumulCasesAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

        # Cumulative cases, control NOT applied at holding over entire simulation
        LogicTempIntVec .= PremVaccStatus.==0 # Control has not been applied.
        LogicTempVec .= PremHasHadInfFlag.==1
        LogicTempIntVec .= LogicTempIntVec.*LogicTempVec #Take product of two conditions. Those that satisfy both will return value 1
        CumulCasesWithoutControlPrem[SimnItr] = sum(LogicTempIntVec)
        LogicTempVec .= LogicTempIntVec.==1
        CumulCasesWithoutControlAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

        # Cumulative cases, control was applied at holding by end of simulation
        LogicTempIntVec .= PremVaccStatus.*PremHasHadInfFlag
        CumulCasesWithControlPrem[SimnItr]  = sum(LogicTempIntVec)
        LogicTempVec .= LogicTempIntVec.==1
        CumulCasesWithControlAnimals[SimnItr,:] = sum(PremLivestockData[LogicTempVec,:],dims=1)

        #-----------------------------------------------------------------------
        #Update node type variables, infectious flag values
        #-----------------------------------------------------------------------
        Infected = findall(IncubationTime .<PremStatus .<=RemovalTime)
        PremInfectiousFlagVec = zeros(Int64,PremNum) #Entry per premises
        PremInfectiousFlagVec[Infected] .= 1  #For those premises infected, amend flag value to 1

        #-----------------------------------------------------------------------
        #Update timing and increment variables
        #-----------------------------------------------------------------------
        t = t + TimeStepVal
        SimnItr = SimnItr + 1

        #-----------------------------------------------------------------------
        #Update EventIndicatorArray & EventTimeArray
        #-----------------------------------------------------------------------

        #Latently infected during timestep
        LogicTempVec .= PremStatus.==TimeStepVal
        EventIndicatorArray[LogicTempVec,2] .= 1
        EventTimeArray[LogicTempVec,2] .= t

        #Infectious during timestep
        LogicTempVec .= PremStatus.==(IncubationTime+TimeStepVal)
        EventIndicatorArray[LogicTempVec,3] .= 1
        EventTimeArray[LogicTempVec,3] .= t

        #Reported during timestep
        LogicTempVec .= PremStatus.==(DetectionTime+TimeStepVal)
        EventIndicatorArray[LogicTempVec,4] .= 1
        EventTimeArray[LogicTempVec,4] .= t

        #Culled during timestep
        LogicTempVec .= CullPremDuringCurrentTimestep.==1
        EventIndicatorArray[LogicTempVec,5] .= 1
        EventTimeArray[LogicTempVec,5] .= t

        #Vaccination events
        if ( (RunControlsFn == CullIPsAndDistanceBasedVaccFn!) ||
            (RunControlsFn == CullIPsAndDistanceBasedCattleVaccFn!) ||
            (RunControlsFn == CullIPsAndRingVaccFn!) ||
            (RunControlsFn == CullIPsAndRingVaccCattleOnlyFn!) ||
            (RunControlsFn == CullIPsAndCumulPremCaseBasedVaccFn!) ||
            (RunControlsFn == CullIPsAndPremPrevalenceBasedVaccFn!) )

            LogicTempVec .= VaccPremDuringCurrentTimestep.==1
            EventIndicatorArray[LogicTempVec,6] .= 1
            EventTimeArray[LogicTempVec,6] .= t

            # Vaccination becomes effective during timestep
            LogicTempVec .= VaccBecomesEffectiveDuringCurrentTimestep.==1
            EventIndicatorArray[LogicTempVec,7] .= 1
            EventTimeArray[LogicTempVec,7] .= t
        end

        #-----------------------------------------------------------------------
        #Exit condition
        #-----------------------------------------------------------------------
        if t > 5
            if (E[SimnItr-3] + I[SimnItr-3] + R2[SimnItr-3]) == 0
                println("E[SimnItr-3]: $(E[SimnItr-3])")
                println("I[SimnItr-3]: $(I[SimnItr-3])")
                println("R2[SimnItr-3]: $(R2[SimnItr-3])")
                IterateFlag = 0
                println("Infection no longer present. Simulation terminated.")
            end
        end

        if mod(SimnItr,20) .== 0
            println("SimnItr: $SimnItr")
        end
    end

    #-------------------------------------------------------------------
    #Write data to files, replicate specific
    #-------------------------------------------------------------------
    writedlm(PremPerDiseaseStatePerTimestepFile, [S E I R R2 Culled VaccPrem CumulCasesPrem CumulCasesWithoutControlPrem CumulCasesWithControlPrem])
    writedlm(CumulativeCulledPerTimestepFile, CullAnimals)
    writedlm(CumulVaccPerTimestepFile, VaccAnimals)
    writedlm(CumulCasesAnimalsPerTimestepFile, [CumulCasesAnimals CumulCasesWithoutControlAnimals CumulCasesWithControlAnimals])

    #-------------------------------------------------------------------
    #Write event indicator arrays to file
    #-------------------------------------------------------------------

    #Find those premises (rows) with non-zero entries for events
    #Column 1 is PremID, so take sum of each row (dims=2), over cols 2:end
    LogicTempIntVec = sum(EventIndicatorArray[:,2:end],dims=2)[:]
    LogicTempVec .= LogicTempIntVec.>0

    #Write only those rows, corresponding to premises where events occurred, to file
    writedlm(EventIndicatorArrayFile, EventIndicatorArray[LogicTempVec,:])
    writedlm(EventTimeArrayFile, EventTimeArray[LogicTempVec,:])

    #-------------------------------------------------------------------
    #Write data to files storing data across all nReps
    #-------------------------------------------------------------------
    writedlm(OutputFileObjs[1], t)
    writedlm(OutputFileObjs[2], [S[SimnItr-1] E[SimnItr-1] I[SimnItr-1] R[SimnItr-1] R2[SimnItr-1] Culled[SimnItr-1] VaccPrem[SimnItr-1] CumulCasesPrem[SimnItr-1] CumulCasesWithoutControlPrem[SimnItr-1] CumulCasesWithControlPrem[SimnItr-1]])
    writedlm(OutputFileObjs[3], CullAnimals[SimnItr-1,:]')
    writedlm(OutputFileObjs[4], VaccAnimals[SimnItr-1,:]')
    writedlm(OutputFileObjs[5], [CumulCasesAnimals[SimnItr-1,:]' CumulCasesWithoutControlAnimals[SimnItr-1,:]' CumulCasesWithControlAnimals[SimnItr-1,:]'])

    #-------------------------------------------------------------------
    #Return outputs (if any)
    #-------------------------------------------------------------------
    return nothing
end
