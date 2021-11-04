#=
Purpose:
Script to compute probability of vaccinated premises not being infected (after moment of vaccination)
Sum over all premises to get expected number of vaccinated premises that would not have been infected

Overview:
 - Function input transmission kernel, premises location and population data
 - Construct premises susceptibility and transmissibility
 - Find premises that underwent vaccination. Construct VaccPresentFlag array.
 - Import premises-level event indicator & time arrays. Construct InfectionPresentFlag & VaccAndInfectionPresentFlag array.
 - Iterate over vaccinated premises.
          -> Get number of days each premises applied infectiousness (uses product of flag arrays)
          -> Compute force of infection for that time period. Account for controls altering transmissibility.
          -> Sum over relevant infectious premises to give cumulative rate.
          -> Compute probability of no infection.
 - Output to file

Julia version: 1.6.3
Date: 3rd November 2021
=#

#--------------------------------------------------------------------------
# LOAD REQUIRED ENVIRONMENT
#--------------------------------------------------------------------------

#Set path to directory this file resides in
cd(dirname(@__FILE__))

using Pkg
Pkg.activate("../../../")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
using DelimitedFiles
using Random
using Distributions

#-------------------------------------------------------------------------------
# IMPORT REQUIRED FUNCTION FILES
#-------------------------------------------------------------------------------
include("../sellke_simn_fns/CalcSuscepTransmissFns.jl")
include("../sellke_simn_fns/DistanceFns.jl")
include("../sellke_simn_fns/SpatialKernels.jl")
include("../sellke_simn_fns/SellkeSimnVarConfigs.jl")

#-------------------------------------------------------------------------------
# DEFINE SUPPORTING FUNCTIONS (Prmeises-level transmission with controls enforced functions)
#-------------------------------------------------------------------------------
"""
    CattleVaccOnlyFn(PremTransmiss_NoControl_SpeciesBreakdown::Array{Float64,2},VaccEff::Float64)

Compute premises transmissibility with vaccination of cattle only.

Inputs:
- `PremTransmiss_NoControl_SpeciesBreakdown::Array{Float64,2}`: Transmission from each livestock type per premises with no controls in place.
- `VaccEff::Float64`: Vaccine effectiveness. Proportion taking value in interval [0,1].

Outputs:
- `PremTransmiss_WithControl::Array{Float64,2}`: Transmission per premises with controls in place.

Location: ComputeVaccButNoInfecStats.jl
"""
function CattleVaccOnlyFn(PremTransmiss_NoControl_SpeciesBreakdown::Array{Float64,2},VaccEff::Float64)

    # Create array equal to dimensions of PremTransmiss_NoControl_SpeciesBreakdown
    PremTransmiss_WithControl_SpeciesBreakdown = 1. *PremTransmiss_NoControl_SpeciesBreakdown

    # Populate array
    PremTransmiss_WithControl_SpeciesBreakdown[:,1] = PremTransmiss_NoControl_SpeciesBreakdown[:,1].*(1-VaccEff) #Cattle impacted by vaccine
    #PremTransmiss_WithControl_SpeciesBreakdown[:,2] = PremTransmiss_NoControl_SpeciesBreakdown[:,2] #Sheep not vaccinated
        # Explicit assignment not needed as initial array equal to PremTransmiss_NoControl_SpeciesBreakdown

    # Sum across columns to get overall premises-level value
    PremTransmiss_WithControl = sum(PremTransmiss_WithControl_SpeciesBreakdown,dims=2)


    return PremTransmiss_WithControl::Array{Float64,2}

end

"""
    AllSpeciesVaccFn(PremTransmiss_NoControl_SpeciesBreakdown::Array{Float64,2},VaccEff::Float64)

Compute premises transmissibility with vaccination of all livestock.

Inputs:
- `PremTransmiss_NoControl_SpeciesBreakdown::Array{Float64,2}`: Transmission from each livestock type per premises with no controls in place.
- `VaccEff::Float64`: Vaccine effectiveness. Proportion taking value in interval [0,1].

Outputs:
- `PremTransmiss_WithControl::Array{Float64,2}`: Transmission per premises with controls in place.

Location: ComputeVaccButNoInfecStats.jl
"""
function AllSpeciesVaccFn(PremTransmiss_NoControl_SpeciesBreakdown::Array{Float64,2},VaccEff::Float64)

    # Create array equal to dimensions of PremTransmiss_NoControl_SpeciesBreakdown
    PremTransmiss_WithControl_SpeciesBreakdown = PremTransmiss_NoControl_SpeciesBreakdown.*(1-VaccEff)

    # Sum across columns to get overall premises-level value
    PremTransmiss_WithControl = sum(PremTransmiss_WithControl_SpeciesBreakdown,dims=2)


    return PremTransmiss_WithControl::Array{Float64,2}

end


#-------------------------------------------------------------------------------
# DEFINE MAIN CALCULATION FUNCTION
#-------------------------------------------------------------------------------
"""
    CalcEstIncorrectVacc(InputFilesPrefix::String,
                            BatchID::String,
                            TotalReplicateNum::Int64,
                            ReplicateOutbreakDurations::Array{Float64,1},
                            PremLivestockData::Array{Int64,2},
                            PremLocs::Array{Float64,2},
                            CoordType::Int64,
                            CalcPremSusceptTransmissFn::Function,
                            TransmissWithControlFn::Function,
                            PerAnimalSuscep::Array{Float64},
                            PerAnimalTransmiss::Array{Float64},
                            SuscepExponent::Array{Float64},
                            TransmissExponent::Array{Float64},
                            KernelFn::Function,
                            delta_t::Float64,
                            VaccEff::Float64,
                            OutputFileName::String)

Calculate the expected number of premises & animals that were vaccinated but would not have been infected if unvaccinated.

Inputs:
- `InputFilesPrefix::String`: Location where premises-level status files reside
- `BatchID::String`: Used as prefix for output files.
- `TotalReplicateNum::Int64`: Replicates performed in total for tested paramter setting
- `ReplicateOutbreakDurations::Array{Float64,1}`: For each replicate, the time period spanned by the simn
- `PremLivestockData::Array{Int64,2}`: Number of each livestock type per premises. Row per premises, column per animal.
- `PremLocs::Array{Float64,2}`: Columns with x/y or lat/long data.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong")
- `CalcPremSusceptTransmissFn::Function`: Use the specified function to calculate premises-level susceptibility and transmissibility
- `TransmissWithControlFn::Function`: Use the specified function to calculate premises-level transmissibility with control enforced.
- `PerAnimalSuscep::Array{Float64}`: Susceptibility scale parameters for each species.
- `PerAnimalTransmiss::Array{Float64}`: Transmissibility scale parameters for each species.
- `SuscepExponent::Array{Float64}`: Susceptibility exponents for each species.
- `TransmissExponent::Array{Float64}`: Transmissibility exponents for each species.
- `KernelFn::Function`: Use the specified kernel form. Defines risk of transmission w.r.t. distance
- `delta_t::Float64`: Timestep between each iteration.
- `VaccEff::Float64`: Vaccine effectiveness. Proportion taking value in interval [0,1]
- `OutputFileName::String`: Location where estimated incorrect vaccination results will be output

Outputs: None \n
Location: ComputeVaccButNoInfecStats.jl
"""
function CalcEstIncorrectVacc(InputFilesPrefix::String,
                                BatchID::String,
                                TotalReplicateNum::Int64,
                                ReplicateOutbreakDurations::Array{Float64,1},
                                PremLivestockData::Array{Int64,2},
                                PremLocs::Array{Float64,2},
                                CoordType::Int64,
                                CalcPremSusceptTransmissFn::Function,
                                TransmissWithControlFn::Function,
                                PerAnimalSuscep::Array{Float64},
                                PerAnimalTransmiss::Array{Float64},
                                SuscepExponent::Array{Float64},
                                TransmissExponent::Array{Float64},
                                KernelFn::Function,
                                delta_t::Float64,
                                VaccEff::Float64,
                                OutputFileName::String)

    # Disaggregate PremLocs
    PremLoc_xVals = PremLocs[:,1]
    PremLoc_yVals = PremLocs[:,2]

    # Compute maximum distances between premises
    MaxDiff_xVals = maximum(PremLoc_xVals)
    MaxDiff_yVals = maximum(PremLoc_yVals)
    MaxDistWithinLandscape = round(Int,sqrt((MaxDiff_xVals*MaxDiff_xVals) + (MaxDiff_yVals*MaxDiff_yVals)))

    # Construct kernel lookup vector
    KernelLookUpVec = KernelFn(MaxDistWithinLandscape)

    # Get number of premises and livestock species under consideration
    # Matches number of rows and columns in PremLivestockData, respectively.
    PremNum = size(PremLivestockData,1)
    LivestockSpeciesNum = size(PremLivestockData,2)
    println("PremNum: $PremNum")
    println("LivestockSpeciesNum: $LivestockSpeciesNum")

    # Initialise storage arrays - estimates for amount of vaccinated entities that would not have been infected
    EstIncorrectVaccPrem = zeros(TotalReplicateNum)
    EstIncorrectVaccAnimals = zeros(TotalReplicateNum,LivestockSpeciesNum)

    # Initialise storage arrays - counts for entities that were both vaccinated AND infected
    NumPremBothVaccAndInfec = zeros(TotalReplicateNum)
    NumAnimalsBothVaccAndInfec = zeros(TotalReplicateNum,LivestockSpeciesNum)

    # Get number of timesteps performed per replicate
    # Add one to acccount for entry at time t=0
    TimestepsPerReplicate = convert(Array{Int64,1},(ReplicateOutbreakDurations./delta_t) .+ 1)

    #-----------------------------------------------------------------------
    ### CONSTRUCT PREMISES SUSCEPTIBILITY AND TRANSMISSIBILITY
    ### THESE VALUES ARE IN THE ABSENCE OF CONTROL
    #-----------------------------------------------------------------------
    PremSuscept_NoControl, PremTransmiss_NoControl,
    PremSuscept_NoControl_SpeciesBreakdown, PremTransmiss_NoControl_SpeciesBreakdown = CalcPremSusceptTransmissFn(PremLivestockData,
                                                                                                                    PerAnimalSuscep, SuscepExponent,
                                                                                                                    PerAnimalTransmiss, TransmissExponent)

    # Iterate over each replicate.
    # Assign expected incorrect vaccination quantities to output arrays/file
    for ReplicateID = 1:TotalReplicateNum

        #PremTransmiss in presence of control
        PremTransmiss_WithControl = TransmissWithControlFn(PremTransmiss_NoControl_SpeciesBreakdown,VaccEff)

        #-----------------------------------------------------------------------
        ### Import premises-level event indicator array
        ### Import premises-level event time array
        #-----------------------------------------------------------------------

        #Array column defs: [PremID  Sus->Latent  Latent->Infectious  Notified  Culled  Vaccinated]
        EventIndicatorArray = readdlm(string(InputFilesPrefix,"EventIndicatorArray_BatchID$(BatchID)_Replicate$(ReplicateID).txt"),Int64)
        EventTimeArray = readdlm(string(InputFilesPrefix,"EventTimeArray_BatchID$(BatchID)_Replicate$(ReplicateID).txt"))


        #Get number of timesteps completed for replicate under consideration
        ReplicateTimeStepTotal = TimestepsPerReplicate[ReplicateID] .- 1 #Subtract one as do not make use with the final state of EventArrays
                                                                        #Final state of event arrays gives status at FinalSimnTime, and would be used up to time FinalSimnTime + delta_t
                                                                        #Therefore, out of time hoirzon undr consideration & not relevant for the force of infection calculaton

        #-----------------------------------------------------------------------
        ### FIND PREMISES THAT UNDERWENT VACCINATION. CONSTRUCT VACCPRESENTFLAG ARRAY
        #-----------------------------------------------------------------------

        #Read final column of EventIndicatorArray
        FinalVaccStatus = EventIndicatorArray[:,end]::Array{Int64,1}

        #Retain rows corresponding to premises that were vaccinated during outbreak
        RetainVaccPremIdxs = FinalVaccStatus.==1
        EventTimeArray_VaccPremRetained = EventTimeArray[RetainVaccPremIdxs,6] #Just keep the values relevant to vaccination time
            # Retain values relevant to time vaccination became effective

        #Get premises IDs that were vaccinated (value equal to 1)
        VaccPrem_AllIDs = EventIndicatorArray[RetainVaccPremIdxs,1]
        TotalVaccPrem = length(VaccPrem_AllIDs)

        #Initialise & construct VaccPresentFlag array
        VaccPresentFlag = zeros(Int64,ReplicateTimeStepTotal,PremNum)
        for VaccPremItr = 1:TotalVaccPrem

            # Get global ID for premises under consideration in this itration of loop
            VaccPremID = VaccPrem_AllIDs[VaccPremItr]

            #Get time at which premises was vaccinated
            #EventTimeArray_VaccPremRetained: Entry per vaccinated premsies.
            VaccTimestepVal = EventTimeArray_VaccPremRetained[VaccPremItr]

            #Get array index corresponding to premises vaccination time (add 1, so time 0 gives index 1 and so on)
            VaccBegins_ArrayIdx = convert(Int64,(VaccTimestepVal/delta_t) + 1)

            #Populate VaccPresentFlag array
            VaccPresentFlag[VaccBegins_ArrayIdx:end,VaccPremID] .= 1
        end

        #-----------------------------------------------------------------------
        ### Find and retain premises with history of being infectious
        #-----------------------------------------------------------------------

        #Retain rows corresponding to premises that were infectious during outbreak
        RetainInfectiousPremIdxs = EventIndicatorArray[:,3].==1
        EventTimeArray_InfectiousPremRetained = EventTimeArray[RetainInfectiousPremIdxs,[2,3,5,6,7]]
            # Extract times for premises becoming
            # [infected, infectious, being culled, undergoing vaccination, vaccination becoming effective]

        #Get premises IDs that were infectious (value equal to 1)
        InfectiousPrem_AllIDs = EventIndicatorArray[RetainInfectiousPremIdxs,1]

        #Get amount of premises infectious during course of outbreak
        #Equates to number of rows of EventIndicatorArray_InfectiousPremRetained
        TotalInfectiousPremNum = size(EventTimeArray_InfectiousPremRetained,1)

        #-----------------------------------------------------------------------
        #INTIALISE & CONSTRUCT INFECTIONPRESENTFLAG & REDUCEDINFECTIONPRESENTFLAG ARRAYS
        #-----------------------------------------------------------------------
        InfectionPresentFlag = zeros(Int64,ReplicateTimeStepTotal,TotalInfectiousPremNum)
        ReducedInfectionPresentFlag = zeros(Int64,ReplicateTimeStepTotal,TotalInfectiousPremNum)
        for InfectiousPremItr = 1:TotalInfectiousPremNum

            # Get time at which premises became infectious & infectious state ended
            # EventTimeArray_InfectiousPremRetained: Row per infectious premsies. Two columns: [TimeBecameInfectious TimeCulled]
            InfectiousStateBegins_TimestepVal = EventTimeArray_InfectiousPremRetained[InfectiousPremItr,2] #Infectiousness begins. Column 2 of EventTimeArray_InfectiousPremRetained
            InfectiousStateEnds_TimestepVal = EventTimeArray_InfectiousPremRetained[InfectiousPremItr,3] #Infectious state over once culled. Column 3 of EventTimeArray_InfectiousPremRetained

            # Get array index corresponding to premises infectious times (add 1, so time 0 gives index 1 and so on)
            InfectiousStateBegins_ArrayIdx = convert(Int64,(InfectiousStateBegins_TimestepVal/delta_t) + 1)

            # Check if premises underwent vaccination. If so, check if successful in reducing tranmission potential
            # (i.e. administered prior to infection being present).
            PremBecameInfectedTime = EventTimeArray_InfectiousPremRetained[InfectiousPremItr,1]
            InfPremVaccTime = EventTimeArray_InfectiousPremRetained[InfectiousPremItr,4]
            InfPremEffecVaccTime = EventTimeArray_InfectiousPremRetained[InfectiousPremItr,5]
            if (InfPremVaccTime == -1) || (InfPremEffecVaccTime == -1) ||
                 (InfPremVaccTime >= PremBecameInfectedTime)
                #Premises never vaccinated or vaccination occurred after premises was infected (so ineffective)

                #Populate InfectionPresentFlag array
                if InfectiousStateEnds_TimestepVal == -1 #Account for situation where premises still infectious at end of simn
                    InfectionPresentFlag[InfectiousStateBegins_ArrayIdx:end,InfectiousPremItr] .= 1
                else #Premises not infectious at end of simulation
                    #Get array index corresponding to premises infectiousness ending (DO NOT ADD 1, want timestep immediately BEFORE culling)
                    InfectiousStateEnds_ArrayIdx = convert(Int64,(InfectiousStateEnds_TimestepVal/delta_t))

                    #Populate InfectionPresentFlag array
                    InfectionPresentFlag[InfectiousStateBegins_ArrayIdx:InfectiousStateEnds_ArrayIdx,InfectiousPremItr] .= 1
                end

            elseif (InfPremEffecVaccTime < PremBecameInfectedTime)
                #Vaccination occurred and was effective prior to infectiousness.

                #Populate InfectionPresentFlag array
                if InfectiousStateEnds_TimestepVal == -1 #Account for situation where premises still infectious at end of simn
                    ReducedInfectionPresentFlag[InfectiousStateBegins_ArrayIdx:end,InfectiousPremItr] .= 1
                else #Premises not infectious at end of simulation
                    #Get array index corresponding to premises infectiousness ending (DO NOT ADD 1, want timestep immediately BEFORE culling)
                    InfectiousStateEnds_ArrayIdx = convert(Int64,(InfectiousStateEnds_TimestepVal/delta_t))

                    #Populate InfectionPresentFlag array
                    ReducedInfectionPresentFlag[InfectiousStateBegins_ArrayIdx:InfectiousStateEnds_ArrayIdx,InfectiousPremItr] .= 1
                end
            else
                # Error check.
                # Should have entered one of the above two loops.
                # Exit programme if did not.
                error("InfPremEffecVaccTime: $InfPremEffecVaccTime. PremBecameInfectedTime: $PremBecameInfectedTime. These values are incompatible.")
            end
        end

        #-----------------------------------------------------------------------
        # ITERATE OVER VACCINATED PREMISES.
        #           -> GET NUMBER OF DAYS EACH PREMISES APPLIED INFECTIOUSNESS.
        #           -> COMPUTE FORCE OF INFECTION FOR THAT TIME PERIOD. ACCOUNT FOR CONTROLS ALTERING TRANSMISSIBILITY.
        #           -> SUM OVER RELEVANT INFECTIOUS PREMISES TO GIVE CUMULATIVE RATE.
        #           -> COMPUTE PROBABILITY OF NO INFECTION.
        #-----------------------------------------------------------------------

        # Get number of premises that were vaccinated and NOT infected
        vacc_and_uninfected_prem_flag = (FinalVaccStatus.*(1 .-EventIndicatorArray[:,2])).==1
            # 1 .-EventIndicatorArray[:,2] returns 0 if premises was ever infected, returns 1 if premises was not infected.
        vacc_and_uninfected_prem_AllIDs = EventIndicatorArray[vacc_and_uninfected_prem_flag,1]
        n_vacc_and_uninfected_prem = length(vacc_and_uninfected_prem_AllIDs)

        #Storage variables used in each iteration
        FlagArray_RefPremVaccAndOtherPremInfectious = zeros(Int64,ReplicateTimeStepTotal,TotalInfectiousPremNum)
        FlagArray_RefPremVaccAndOtherPremReducedInfectious = zeros(Int64,ReplicateTimeStepTotal,TotalInfectiousPremNum)

        InfectiousTimestepsPerPremAgainstVaccRefPrem = zeros(Int64,TotalInfectiousPremNum)
        ReducedInfectiousTimestepsPerPremAgainstVaccRefPrem = zeros(Int64,TotalInfectiousPremNum)

        #Output variables
        VaccPrem_NoInfecProb = zeros(n_vacc_and_uninfected_prem)
        VaccAnimals_NoInfecProb = zeros(n_vacc_and_uninfected_prem,LivestockSpeciesNum)

        #Enter iteration over vaccinated premises
        for VaccAndUninfecPremIdx = 1:n_vacc_and_uninfected_prem

            #Get premises ID for vaccinated premises within this loop iteration
            VaccAndUninfecPremID = vacc_and_uninfected_prem_AllIDs[VaccAndUninfecPremIdx]

            #For this reference premises, get vaccination status over time
            RefPrem_VaccPresentFlag = VaccPresentFlag[:,VaccAndUninfecPremID]

            #Multiply vaccination present in reference premises by infectiousness state of all other premises
            FlagArray_RefPremVaccAndOtherPremInfectious .= InfectionPresentFlag.*RefPrem_VaccPresentFlag
            FlagArray_RefPremVaccAndOtherPremReducedInfectious .= ReducedInfectionPresentFlag.*RefPrem_VaccPresentFlag

            #For other premises, get number of infectiousness days applied to (in vaccinated state) reference premises
            #Sum over each columns of FlagArray_RefPremVaccAndOtherPremInfectious (for unvaccinated infectious premises)
            #               & InfectiousTimestepsPerVaccInfPremAgainstVaccRefPrem (for vaccianted infectious premises, that will have modified transmission)
            InfectiousTimestepsPerPremAgainstVaccRefPrem .= sum(FlagArray_RefPremVaccAndOtherPremInfectious,dims=1)[:]
            ReducedInfectiousTimestepsPerPremAgainstVaccRefPrem .= sum(FlagArray_RefPremVaccAndOtherPremReducedInfectious,dims=1)[:]

            #Compute cumulative force of infection against VaccAndUninfecPremID, after vaccination for remainder of simulation
            CumulFOI = 0. #Initialise value
            for InfectiousPremItr = 1:TotalInfectiousPremNum #Iterate over infectious premises.

                #Check infectious premises is not also the reference premises!
                #Check number of combined timesteps (overlap between infectious window and reference prem being vaccinated) exceeds 0.
                #Otherwise, go to next InfectiousPremItr
                OverlapInfectiousVaccTimesteps = InfectiousTimestepsPerPremAgainstVaccRefPrem[InfectiousPremItr]
                OverlapReducedInfectiousVaccTimesteps = ReducedInfectiousTimestepsPerPremAgainstVaccRefPrem[InfectiousPremItr]

                InfectiousPremID = InfectiousPrem_AllIDs[InfectiousPremItr]
                if (OverlapInfectiousVaccTimesteps+OverlapReducedInfectiousVaccTimesteps > 0) && (InfectiousPremID != VaccAndUninfecPremID)

                    #Calculate distance between the infectious premises and reference premises
                    if CoordType == 1 #Cartesian (metres unit)
                        xDiff = PremLoc_xVals[InfectiousPremID] - PremLoc_xVals[VaccAndUninfecPremID]
                        yDiff = PremLoc_yVals[InfectiousPremID] - PremLoc_yVals[VaccAndUninfecPremID]
                        D = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                    elseif CoordType == 2 #Cartesian (km unit)
                        xDiff = PremLoc_xVals[InfectiousPremID] - PremLoc_xVals[VaccAndUninfecPremID]
                        yDiff = PremLoc_yVals[InfectiousPremID] - PremLoc_yVals[VaccAndUninfecPremID]
                        D = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                    elseif CoordType == 3 #LatLong
                        xDiff = PremLoc_xVals[InfectiousPremID] - PremLoc_xVals[VaccAndUninfecPremID]
                        yDiff = PremLoc_yVals[InfectiousPremID] - PremLoc_yVals[VaccAndUninfecPremID]
                        D = sqrt((xDiff*xDiff) + (yDiff*yDiff))
                    end

                    #Amend distance to an integer in metres
                    distIdx = ReturnDistIdxForKernel(D)

                    #Get total elapsed time. Timesteps multiplied by timestep length.
                    ElapsedTime_FullTransmissFOI = OverlapInfectiousVaccTimesteps*delta_t
                    ElapsedTime_ReducedTransmissFOI = OverlapReducedInfectiousVaccTimesteps*delta_t #When infectious premises has also been succesfully vaccinated

                    #Calculate force of infection for infectious site against vaccinated premises
                    CurrentInfEventFOI = PremSuscept_NoControl[VaccAndUninfecPremID]*KernelLookUpVec[distIdx]*
                                            ((PremTransmiss_NoControl[InfectiousPremID]*ElapsedTime_FullTransmissFOI) +
                                                (PremTransmiss_WithControl[InfectiousPremID]*ElapsedTime_ReducedTransmissFOI))

                    #Update CumulFOI
                    CumulFOI = CumulFOI + CurrentInfEventFOI
                end
            end

            # Get probability of no infection for premises where control applied and that
            # were NOT infected during the simulation
            # (prob. no event in interval with rate lambda is exp(-lambda))
            NoInfecProb = exp(-CumulFOI)

            #Check probability value is valid
            if (NoInfecProb < 0) || (NoInfecProb > 1)
                error("NoInfecProb is invalid. Has value $NoInfecProb. Should be in interval [0,1].")
            end

            #Assign no infection probability to file output variables
            VaccPrem_NoInfecProb[VaccAndUninfecPremIdx] =  NoInfecProb
            VaccAnimals_NoInfecProb[VaccAndUninfecPremIdx,:] = PremLivestockData[VaccAndUninfecPremID,:].*NoInfecProb
        end

        #-----------------------------------------------------------------------
        ### ASSIGN ESTIMATED UNNEEDED VACCINATION FOR REPLICATE TO STORAGE VECTOR
        #-----------------------------------------------------------------------
        EstIncorrectVaccPrem[ReplicateID] = sum(VaccPrem_NoInfecProb)
        EstIncorrectVaccAnimals[ReplicateID,:] = sum(VaccAnimals_NoInfecProb,dims=1) #Sum columns to get value for each livestock type

        #-----------------------------------------------------------------------
        ### COMPUTE NUMBER OF PREMISES/ANIMALS BOTH VACCINATED AND INFECTED
        #-----------------------------------------------------------------------
        #Product of ever being infected vector and ever being vaccinated vector. Entries of 1 if both satisfied
        PremBothVaccAndInfecFlag = EventIndicatorArray[:,2].*FinalVaccStatus
        NumPremBothVaccAndInfec[ReplicateID] = sum(PremBothVaccAndInfecFlag)

        #Get livestock data for premises where events occurred
        #Then compute the number of animals both vaccinated and infected
        PremsWithEventsOccurring_AllIDs = EventIndicatorArray[:,1] #Obtain IDs for premises where events occurred
        PremsWithEventsOccurring_LivestockData = PremLivestockData[PremsWithEventsOccurring_AllIDs,:]
        NumAnimalsBothVaccAndInfec[ReplicateID,:] = sum(PremsWithEventsOccurring_LivestockData[PremBothVaccAndInfecFlag.==1,:],dims=1) #Sum columns to get value for each livestock type

        #Output progress to screen
        println("ReplicateID $ReplicateID complete.")
    end

    #---------------------------------------------------------------------------
    ### OUTPUT TO FILE
    #---------------------------------------------------------------------------

    #Estimates for amount of vaccinated entities that would not have been infected
    writedlm(string(OutputFileName,"EstIncorrectVaccPrem_BatchID$(BatchID).txt"),EstIncorrectVaccPrem)
    writedlm(string(OutputFileName,"EstIncorrectVaccAnimals_BatchID$(BatchID).txt"),EstIncorrectVaccAnimals)

    #Counts for entities that were both vaccinated AND infected
    writedlm(string(OutputFileName,"NumPremBothVaccAndInfec_BatchID$(BatchID).txt"),NumPremBothVaccAndInfec)
    writedlm(string(OutputFileName,"NumAnimalsBothVaccAndInfec_BatchID$(BatchID).txt"),NumAnimalsBothVaccAndInfec)

    #---------------------------------------------------------------------------
    ### RETURN VARIABLES FROM FUNCTION
    #---------------------------------------------------------------------------
    return nothing
end

#--------------------------------------------------------------------------
# SET VARIABLES FROM ARGS
#--------------------------------------------------------------------------
args = ARGS

# If running locally from REPL, will not have any input ARGS
# Can set values for input parameters to main run function here
# args/ARGS list
# args[1] job_ID
# args[2] job_ID_foundation_value: Set value that job_ID will be incremented upon
# args[3] BatchID_offset: Value that BatchID is offset by
# args[4] ConfigFn: Location and pathogen configuration
# args[5] nReps: Number of replicates requested
if length(ARGS)==0
    args = [ "11","0","0","CumbriaSellkeFMDconfig","10"]
end

# To run from command line, example:
# julia ComputeVaccButNoInfecStats.jl 1 0 CumbriaFMDconfig 10

# Set identifier for job
const job_ID = parse(Int64, args[1])

# Set value that job_ID will be incremented upon
const job_ID_foundation_value = parse(Int64, args[2])

# Set value that BatchID will be offset by
const BatchID_offset = parse(Int64, args[3])

# # Specify configuration function
# # Relevant files in "../GridSimnFns/SimnVarConfigs.jl"
const ConfigFn = args[4] #Make Symbol a callable function

# State the number of replicates that were run
const TotalReplicateNum = parse(Int64, args[5])

# Set range of BatchIDs to be processed
# const BatchID_start_idx = BatchID_offset + (job_ID-1)*100 + 1
# const BatchID_end_idx = BatchID_offset + (job_ID)*100
const BatchID_start_idx = BatchID_offset + job_ID_foundation_value + job_ID
const BatchID_end_idx = BatchID_offset + job_ID_foundation_value + job_ID

# Specify BatchID to be checked
const BatchID_vec = collect(BatchID_start_idx:1:BatchID_end_idx)

#--------------------------------------------------------------------------
# CALL FUNCTION WITH SPECIFIED VARIABLE VALUES
#--------------------------------------------------------------------------

# Set input file name prefix
if (ConfigFn == "CumbriaSellkeFMDconfig") ||
        (ConfigFn == "CumbriaSellkeFMDconfig_alternate_transmiss_params")
    InputFilesPrefix = "../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_EventArrays/"
elseif (ConfigFn == "DevonSellkeFMDconfig") ||
        (ConfigFn == "DevonSellkeFMDconfig_alternate_transmiss_params")
    InputFilesPrefix = "../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_EventArrays/"
elseif (ConfigFn == "TestSellkeSimnConfig")
    InputFilesPrefix = "../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs_EventArrays/"
else
    error("Invalid ConfigFn input.")
end

for batch_itr = 1:length(BatchID_vec)

    # BatchID - Used as prefix for output files
    BatchID = string(BatchID_vec[batch_itr])

    if (ConfigFn == "CumbriaSellkeFMDconfig") ||
            (ConfigFn == "CumbriaSellkeFMDconfig_alternate_transmiss_params")
        # ReplicateOutbreakDurations - For each replicate, the time period spanned by the simn
        ReplicateOutbreakDurations = readdlm("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/OutbreakDuration_BatchID$(BatchID).txt")
        ReplicateOutbreakDurations = ReplicateOutbreakDurations[:] #Pass to function as 1D vector

        ### Load Cumbria livestock and location data ###
        CountyRawData = readdlm("../../../data/dummy_data/cumbria_synthetic_data.txt")
        # CountyRawData = readdlm("../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID08.txt")
    elseif (ConfigFn == "DevonSellkeFMDconfig") ||
            (ConfigFn == "DevonSellkeFMDconfig_alternate_transmiss_params")
        # ReplicateOutbreakDurations - For each replicate, the time period spanned by the simn
        ReplicateOutbreakDurations = readdlm("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/OutbreakDuration_BatchID$(BatchID).txt")
        ReplicateOutbreakDurations = ReplicateOutbreakDurations[:] #Pass to function as 1D vector

        ### Load Devon livestock and location data ###
        CountyRawData = readdlm("../../../data/dummy_data/devon_synthetic_data.txt")
        # CountyRawData = readdlm("../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID10.txt")
    elseif (ConfigFn == "TestSellkeSimnConfig")
        # ReplicateOutbreakDurations - For each replicate, the time period spanned by the simn
        ReplicateOutbreakDurations = readdlm("../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs_Aggregated/OutbreakDuration_BatchID$(BatchID).txt")
        ReplicateOutbreakDurations = ReplicateOutbreakDurations[:] #Pass to function as 1D vector

        ### Load dummy livestock and location data ###
        CountyRawData = readdlm("../../../data/dummy_data/livestock_dummy_data.txt")
    else
        error("Invalid ConfigFn. Could not load CountyRawData.")
    end

    # Assign cattle herd and sheep flock size to variable
    # Retain only Cattle and Sheep (columns 5 and 6)
    # Cast as integer type array
    PremLivestockDataAll = convert(Array{Int64,2},abs.(CountyRawData[:,[5,6]]))
        #Negative values in data designate user-input value. Take absolute value

    # Retain those premises that have cattle and/or sheep present
    TotalCattleSheepOnPrem = sum(PremLivestockDataAll,dims=2)
    RetainPremFlag = TotalCattleSheepOnPrem[:].>0 #Use [:] to ensure RetainPremFlag is type BitArray{1}, not BitArray{2}!
    PremLivestockData = PremLivestockDataAll[RetainPremFlag,:]

    #Assign location info to variable
    #Lat & Long co-ordinates correspond to columns 3 and 4 of input data
    PremLocs = convert(Array{Float64,2},abs.(CountyRawData[RetainPremFlag,3:4])) #Negative values in data designate user-input value. Take absolute value
    ################################################################################

    # # If needed, remove premises that lie outside the specified bounding box
    # if ( (ConfigFn == "DevonSellkeFMDconfig") ||
    #      (ConfigFn == "DevonSellkeFMDconfig_alternate_transmiss_params") )
    #
    #      # Devon region bounding box
    #      BoundingBoxLeft = 212000.
    #      BoundingBoxRight = 340000.
    #      BoundingBoxBottom = 34960.
    #      BoundingBoxTop = 151210.
    #
    #      # Check if locations lie inside bounding box
    #      valid_prem_loc_flag = (PremLocs[:,1] .>= BoundingBoxLeft) .&
    #                              (PremLocs[:,1] .<= BoundingBoxRight) .&
    #                              (PremLocs[:,2] .>= BoundingBoxBottom) .&
    #                              (PremLocs[:,2] .<= BoundingBoxTop)
    #
    #      # Update location and livestock arrays
    #      PremLocs  = PremLocs[valid_prem_loc_flag,:]
    #      PremLivestockData = PremLivestockData[valid_prem_loc_flag,:]
    # end


    ################################################################################

    # CoordType - Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong")
    CoordType = 1

    ### Suscept & Transmiss related variables ###
    CalcPremSusceptTransmissFn = CalcPremSuscepTransmiss_FMDlike

    #TransmissWithControlFn = CattleVaccOnlyFn
    TransmissWithControlFn = AllSpeciesVaccFn

    if (ConfigFn == "CumbriaSellkeFMDconfig")
        # Cumbria parameter set
        PerAnimalSuscep = [5.7 1] #From 2008 FMD paper for Cumbria
        PerAnimalTransmiss = [8.2e-4 8.3e-4] #From 2008 FMD paper for Cumbria
        SuscepExponent = [0.41 0.2] #From 2008 FMD paper for Cumbria
        TransmissExponent = [0.42 0.49] #From 2008 FMD paper for Cumbria
    elseif (ConfigFn == "CumbriaSellkeFMDconfig_alternate_transmiss_params")
        # Cumbria parameter set, reduced per animal tranmissibility and elongated infectious period
        PerAnimalSuscep = [5.7 1] #From 2008 FMD paper for Cumbria
        PerAnimalTransmiss = [8.2e-4 8.3e-4].*0.25 #From 2008 FMD paper for Cumbria
        SuscepExponent = [0.41 0.2] #From 2008 FMD paper for Cumbria
        TransmissExponent = [0.42 0.49] #From 2008 FMD paper for Cumbria
    elseif (ConfigFn == "DevonSellkeFMDconfig")
        PerAnimalSuscep = [5.7 1]
        PerAnimalTransmiss = [8.2e-4 8.3e-4]
        SuscepExponent = [0.41 0.2]
        TransmissExponent = [0.42 0.49]
    elseif (ConfigFn == "DevonSellkeFMDconfig_alternate_transmiss_params")
        PerAnimalSuscep = [5.7 1]
        PerAnimalTransmiss = [8.2e-4 8.3e-4].*0.25
        SuscepExponent = [0.41 0.2]
        TransmissExponent = [0.42 0.49]
    elseif (ConfigFn == "TestSellkeSimnConfig")
        PerAnimalSuscep = [5.7 1]
        PerAnimalTransmiss = [8.2e-4 8.3e-4].*0.25
        SuscepExponent = [0.41 0.2]
        TransmissExponent = [0.42 0.49]
    end
    ################################################################################

    # KernelFn - (function name) Use the specified kernel form. Defines risk of transmission w.r.t. distance
    # Kernel options housed in "CommonSimnFns/SpatialKernels.jl"
    KernelFn = Construct_USDOS2_kernel

    #delta_t - Timestep between each iteration.
    delta_t = 1.

    #Specify vaccine efficacy used in the simulations
    VaccEff = 1.

    # OutputFileName - Location where estimated incorrect vaccination results will be output
    if (ConfigFn == "CumbriaSellkeFMDconfig") ||
            (ConfigFn == "CumbriaSellkeFMDconfig_alternate_transmiss_params")
        # Cumbria runs, Sellke construction
        OutputFileName = "../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/"
    elseif (ConfigFn == "DevonSellkeFMDconfig") ||
            (ConfigFn == "DevonSellkeFMDconfig_alternate_transmiss_params")
        # Devon runs, Sellke construction
        OutputFileName = "../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/"
    elseif (ConfigFn == "TestSellkeSimnConfig")
        # Dummy model, Sellke construction
        OutputFileName = "../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs_Aggregated/"
    end

    #Value error checks
    if (CoordType != 1) && (CoordType != 2) && (CoordType != 3)
        error("CoordType has value $VaccEff. Needs to take value 1, 2 or 3.")
    end

    if (VaccEff < 0) || (VaccEff > 1) #Error check
        error("VaccEff has value $VaccEff. Needs to be in interval [0,1].")
    end

    #Call function
    CalcEstIncorrectVacc(InputFilesPrefix,
                                    BatchID,
                                    TotalReplicateNum,
                                    ReplicateOutbreakDurations,
                                    PremLivestockData,
                                    PremLocs,
                                    CoordType,
                                    CalcPremSusceptTransmissFn,
                                    TransmissWithControlFn,
                                    PerAnimalSuscep,
                                    PerAnimalTransmiss,
                                    SuscepExponent,
                                    TransmissExponent,
                                    KernelFn,
                                    delta_t,
                                    VaccEff,
                                    OutputFileName)
    println("Iteration $batch_itr complete")
end
