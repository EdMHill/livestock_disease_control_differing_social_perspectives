#=
Purpose:
Run an individual-based livestock disease epidemic model using the Sellke construction

Julia version: 1.6.3
Date: 3rd November 2021
=#

#--------------------------------------------------------------------------
# LOAD ENVIRONMENT
#--------------------------------------------------------------------------

#Set path to directory this file resides in
cd(dirname(@__FILE__))

using Pkg
Pkg.activate("../../../")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
using DelimitedFiles
using Distributions
using Statistics
using StatsBase
using Random
using JLD2
using Shapefile
using Profile

#-------------------------------------------------------------------------------
# SET UP COMPOSITE TYPES
#-------------------------------------------------------------------------------

#Landscape information
struct LandscapeData
    get_n_nodes::Int64  #Number of nodes in the landscaspe
    longest_side::Float64  #Length of longest dimension of box containing all premises
    get_trans::Float64
    get_susc::Float64
    coordType::Int64 #Type of coordinate system in use
end

#-------------------------------------------------------------------------------
# IMPORT REQUIRED FUNCTION FILES
#-------------------------------------------------------------------------------
include("../sellke_simn_fns/CommonFunctions.jl")
include("../sellke_simn_fns/Grid.jl") # For SeedInfection fn & ReturnDistIdxForKernel
include("../sellke_simn_fns/CalcSuscepTransmissFns.jl")
include("../sellke_simn_fns/SpatialKernels.jl")
include("../sellke_simn_fns/ControlFns.jl")
include("../sellke_simn_fns/DistanceFns.jl")
include("../sellke_simn_fns/inoculation_time_fns.jl")
include("../sellke_simn_fns/SellkeSimnVarConfigs.jl")
include("../sellke_simn_fns/SellkeLocalSpreadFns.jl")
include("../sellke_simn_fns/RunSellkeOutbreakReplicateFns.jl")

#-------------------------------------------------------------------------------
# FUNCTIONS FOR RUNNING SPATIAL MODEL
#-------------------------------------------------------------------------------

# Wrapper to run outbreak model that uses the Sellke construction
# Main chunks: (i) Disaggregate input tuples
#              (ii) Construct kernel lookup vector
#              (iii) Get node transmissibility & susceptibility
#              (iv) If not calculating distances on the fly, construct distances array,
#                     apply kernel and calculate FOI array
#               (v) Run main iteration. Calls function to run a single replicate
"""
    SpatialLivestockSellkeInfectionModel(BatchID::String,
                                        PremLocData::Array{Any,1},
                                        PremLivestockData::Array{Int64},
                                        TimeParams::Array{Float64,1},
                                        nReps::Int64,
                                        CalcPremSuscepTransmissFn::Function,
                                        KernelFn::Function,
                                        GridOptimParams, #Type of each entry in tuple is cast later in script
                                        EpiParamVals::Array{Float64,1},
                                        ControlParamVals::Array{Any,1},
                                        LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                        InitialInfInfo, #Type of each entry in tuple is cast within SeedInfection function (file Grid.jl)
                                        OutputLevels, #Either Int64, value 0, or an Array{Int64,1}.
                                        PerAnimalSuscep::Array{Float64},
                                        PerAnimalTransmiss::Array{Float64},
                                        SuscepExponent::Array{Float64},
                                        TransmissExponent::Array{Float64},
                                        IterateOutbreakFn::Function,
                                        RunControlsFn::Function,
                                        RunSimulationReplicateFn::Function,
                                        SaveGridConfigFileNames::Array{String,1},
                                        ReplicateFilePrefix::String,
                                        OutputFileObjs::Array{IOStream,1})

Wrapper to run outbreak model on selected configuration using the Sellke costruction.

Inputs:
- `BatchID::String`: Used as prefix for output files.
- `PremLocData::Array{Any,1}`: First entry: Co-ordinate type; Second entry: Columns with x/y or lat/long data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
- `PremLivestockData::Array{Int64}`: Number of each livestock type per premises. Row per premises, column per animal.
- `TimeParams::Array{Float64,1}`: Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over.
- `nReps::Int64`: Number of replicates to run.
- `CalcPremSuscepTransmissFn::Function`: Use the specified function to calculate premises-level susceptibility and transmissibility
- `KernelFn::Function`: Use the specified kernel form. Defines risk of transmission w.r.t. distance.
- `GridOptimParams`: Entry one: (flag) Specify whether grid optimisation method should be used. (1 = fixed, 2 = dynamic). Entry two: (integer) Corresponds to either: when first entry = 1, number of cells along one dimension of the grid (total number will be this^2). when first entry = 2, nMax number of nodes in one cell when using dynamic cell size. Type of each entry in tuple is cast later in script.
- `EpiParamVals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)]
- `ControlParamVals::Array{Any,1}`: Variables related to implementing control measures
- `LinkedPremToControlIdxs::Array{Array{Int64,1}}`: For each premises, those linked premises that would also undergo control
- `InitialInfInfo`: Type of each entry in tuple is cast within SeedInfection function (file Grid.jl).  [SeedMethod,NumOfNodes/NodeIDs]
                             Seed method. 1 = random,
                                           2 = single specific node id,
                                           3 = set of specific node ids,
                                           4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
                             Number of nodes to seed each replicate if SeedMethod = 1; or if SeedMethod = 2/3, seed these specific node every replicate.
                             Node ids to seed. One id/row, number of lines must be == number of replicates.
- `OutputLevels`: Will output replicate results when the cumulative number of infected reaches these levels. Will also exit after the last output. Set to 0 if want to output every timestep.
- `PerAnimalSuscep::Array{Float64}`: Susceptibility scale parameters for each species.
- `PerAnimalTransmiss::Array{Float64}`: Transmissibility scale parameters for each species.
- `SuscepExponent::Array{Float64}`: Susceptibility exponents for each species.
- `TransmissExponent::Array{Float64}`: Transmissibility exponents for each species.
- `IterateOutbreakFn::Function`: Function to perform disease transitions and control implementation per timestep
- `RunControlsFn::Function`: Enact controls as specified within the given function.
- `RunSimulationReplicateFn::Function`: Function to perform single outbreak replicate
- `SaveGridConfigFileNames::Array{String,1}`: Names for three files. (i) Save array defining boundary limits of each cell within obtained grid configuration
                                  Array columns correspond to [xMin, xMax, yMin, yMax];
                              (ii) Grid ID each premises resides in;
                              (iii) For each grid ID, total number of premises in the cell.
- `ReplicateFilePrefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `OutputFileObjs::Array{IOStream,1}`: Filename identifiers written to by all replicates.
- `precalc_dist_and_kernel_array_flag::Bool`: Indicator that if true, construct array of premises-to-premises kernel quantities before running main outbreak function.

Outputs: None \n
Location: run\\_generic\\_livestock\\_disease\\_spatial\\_model.jl
"""
function SpatialLivestockSellkeInfectionModel(BatchID::String,
                                            PremLocData::Array{Any,1},
                                            PremLivestockData::Array{Int64},
                                            TimeParams::Array{Float64,1},
                                            nReps::Int64,
                                            CalcPremSuscepTransmissFn::Function,
                                            KernelFn::Function,
                                            EpiParamVals::Array{Float64,1},
                                            ControlParamVals::Array{Any,1},
                                            LinkedPremToControlIdxs::Array{Array{Int64,1}},
                                            InitialInfInfo, #Type of each entry in tuple is cast within SeedInfection function (file Grid.jl)
                                            OutputLevels, #Either Int64, value 0, or an Array{Int64,1}.
                                            PerAnimalSuscep::Array{Float64},
                                            PerAnimalTransmiss::Array{Float64},
                                            SuscepExponent::Array{Float64},
                                            TransmissExponent::Array{Float64},
                                            IterateOutbreakFn::Function,
                                            RunControlsFn::Function,
                                            RunSimulationReplicateFn::Function,
                                            ReplicateFilePrefix::String,
                                            OutputFileObjs::Array{IOStream,1},
                                            precalc_dist_and_kernel_array_flag::Bool)

    #---------------------------------------------------------------------------
    # (I) DISAGGREAGTE INPUT TUPLES
    #---------------------------------------------------------------------------

    #Time related variables
    TimeStepVal = TimeParams[1]::Float64
    MaxTime = TimeParams[2]::Float64

    #Epidemioloigcal parameters
    IncubationTime = EpiParamVals[1]::Float64
    DetectionTime = EpiParamVals[2]::Float64
    RemovalTime = EpiParamVals[3]::Float64

    #Co-ordinate type variable
    CoordType = PremLocData[1]::Int64
    if CoordType != 1 && CoordType != 2 && CoordType != 3     #Check value, throw error if invalid
        error("CoordType has value $CoordType, but CoordType must take value 1 ('Cartesian (metres)'), 2 ('Cartesian (km)') or 3 ('LatLong').")
    end
    # if CoordType != "Cartesian" && CoordType != "LatLong"     #Check value, throw error if invalid
    #     error("CoordType has value $CoordType, but CoordType must take value 'Cartesian' or 'LatLong'.")
    # end

    #Node location attributes
    PremNum = size(PremLocData[2],1)::Int64 #Rows of premises location dataset, stored in entry 2 of PremLocData tuple
    PremLoc_AllVals = PremLocData[2]::Array{Float64,2}
    PremLoc_xVals = PremLocData[2][:,1]::Array{Float64,1} #Extract positional data from PremLocData tuple
    PremLoc_yVals = PremLocData[2][:,2]::Array{Float64,1}
    #If co-ordinate type one, columns are x/y data
    #If co-ordinate type two, columns are lat/long data

    #---------------------------------------------------------------------------
    # (II) CONSTRUCT KERNEL LOOKUP VECTOR
    #---------------------------------------------------------------------------
    #Disaggregate landscape bounding box variables.
    BoundingBoxVar = PremLocData[3]::Array{Float64,1}   #Vector with entries [Min_x,Max_x,Min_y,Max_y]

    #Call function to get landscape dimensions
    BoundingBoxWidth, BoundingBoxHeight, LongestLandscapeEdge, MaxDistWithinLandscape =
                    GetLandscapeSize(BoundingBoxVar,
                                        PremLoc_xVals,
                                        PremLoc_yVals,
                                        CoordType)
    # Get the kernel lookup vector
    KernelLookUpVec = KernelFn(MaxDistWithinLandscape)

    #---------------------------------------------------------------------------
    # (III) GET NODE TRANSMISSIBILITY & SUSCEPTIBILITY
    # ALSO GET SPECIES LEVEL TRANSMISSIBILITY & SUSCEPTIBILITY (PER NODE)
    #---------------------------------------------------------------------------

    # Call function to calculate premises susceptibilities and transmissability
    # that will be used to reinitialise variables in each replicate
    PremSuscept_start_val, PremTransmiss_start_val,
    PremSuscept_SpeciesBreakdown_start_val, PremTransmiss_SpeciesBreakdown_start_val = CalcPremSuscepTransmissFn(PremLivestockData,
                                                                                                    PerAnimalSuscep, SuscepExponent,
                                                                                                    PerAnimalTransmiss, TransmissExponent)

    #---------------------------------------------------------------------------
    # (IV) IF NOT CALCULATING DISTANCES ON THE FLY, CONSTRUCT DISTANCES ARRAY & APPLY KERNEL
    #---------------------------------------------------------------------------
    if precalc_dist_and_kernel_array_flag == true

        # Initialise kernel_val_array
        kernel_val_array = zeros(Float64,PremNum,PremNum)

        # Iterate over each premises combination & populate kernel_val_array
        for prem_ii_ID = 1:PremNum
            # Get location of premises ii
            prem_ii_xLoc = PremLoc_xVals[prem_ii_ID]
            prem_ii_yLoc = PremLoc_yVals[prem_ii_ID]

            for prem_jj_ID = prem_ii_ID:PremNum # Matrix is symmetric
                # Get location of premises jj
                prem_jj_xLoc = PremLoc_xVals[prem_jj_ID]
                prem_jj_yLoc = PremLoc_yVals[prem_jj_ID]

                #Calculate distance between the two points
                if CoordType == 1 #Cartesian co-ords (metres)
                    d =  eucl_distance(prem_ii_xLoc,
                                        prem_ii_yLoc,
                                        prem_jj_xLoc,
                                        prem_jj_yLoc)
                elseif CoordType == 2 #Cartesian co-ords (metres)
                    d = eucl_distance_ConvertToMetres(prem_ii_xLoc,
                                                        prem_ii_yLoc,
                                                        prem_jj_xLoc,
                                                        prem_jj_yLoc)
                elseif CoordType == 3 #Lat/Long co-ords
                    d = GreatCircleDistance(prem_ii_yLoc, prem_ii_xLoc,  #lat1, lon1
                                            prem_jj_yLoc, prem_jj_xLoc) #lat2, lon2
                end

                # Get distance between the two premises
                distIdx = ReturnDistIdxForKernel(d)

                # Apply the kernel to that distance
                kernel_val_array[prem_ii_ID,prem_jj_ID] = KernelLookUpVec[distIdx]
                kernel_val_array[prem_jj_ID,prem_ii_ID] = KernelLookUpVec[distIdx]
            end
        end
    elseif precalc_dist_and_kernel_array_flag == false
        kernel_val_array = Array{Float64,2}(undef,0,0)
    else
        error("precalc_dist_and_kernel_array_flag has value $(precalc_dist_and_kernel_array_flag). Must take value true or false.")
    end

    #---------------------------------------------------------------------------
    # (V) RUN MAIN ITERATION. CALLS FUNCTION TO RUN A SINGLE REPLICATE
    #---------------------------------------------------------------------------
    SpeciesNum = size(PremLivestockData,2) #Get number of species/livestock types in use
    for ItrIdx = 1:nReps     #Run nReps total replicates

        #-----------------------------------------------------------------------
        # SET UP RANDOM NUMBER GENERATOR AT START OF REPLICATE
        #-----------------------------------------------------------------------
        rng = MersenneTwister(RNGseed+ItrIdx)

        #-----------------------------------------------------------------------
        # REINITIALISE STATUS VARIABLES
        #-----------------------------------------------------------------------

        #Initialise PremStatus, PremHasHadInfFlag and PremVaccStatus
        PremStatus = zeros(PremNum)
        PremHasHadInfFlag = zeros(Int64, PremNum)
        PremVaccStatus = zeros(Int64, PremNum)
        SpeciesGroupVaccStatusByPrem = zeros(Int64, PremNum, SpeciesNum)

        # Reset premises-level and livestock-level susceptiblity and transmissibility
        # to unmodified values (when no controls applied)
        PremSuscept = copy(PremSuscept_start_val)
        PremTransmiss = copy(PremTransmiss_start_val)
        PremSuscept_SpeciesBreakdown = copy(PremSuscept_SpeciesBreakdown_start_val)
        PremTransmiss_SpeciesBreakdown = copy(PremTransmiss_SpeciesBreakdown_start_val)

        #-----------------------------------------------------------------------
        # SEED INFECTION
        #-----------------------------------------------------------------------

        #Initialise disease status indicator array and event time array. Events depends on control function in use
        if ( (RunControlsFn == CullIPsAndDistanceBasedVaccFn!) ||  # Vaccination controls in use
                (RunControlsFn == CullIPsAndDistanceBasedCattleVaccFn!) ||
                (RunControlsFn == CullIPsAndRingVaccFn!) ||
                (RunControlsFn == CullIPsAndRingVaccCattleOnlyFn!) )
            StatesToTrack = length(EpiParamVals) + 3
                # Movements between each disease state from latently infected onward, length(EpiParamVals),
                # & becoming infected event (+1)
                # & implementation of vaccination control (+2, vaccination administered time & vaccination becomes effective time)
        else #Culling controls only
            StatesToTrack = length(EpiParamVals) + 1
            # Movements between each disease state from latently infected onward, length(EpiParamVals),
            # & becoming infected event (+1)
        end
        EventIndicatorArray = [1:1:PremNum zeros(Int64,PremNum,StatesToTrack)] #First column is PremID.
        EventTimeArray = [1:1:PremNum  -1 .*ones(PremNum,StatesToTrack)]

        #Call function to seed infection
        SeedPremIDs =  SeedInfection(InitialInfInfo,
                                            PremNum::Int64,
                                            PremLoc_AllVals,
                                            CoordType,
                                            rng)
        println("SeedPremIDs: $SeedPremIDs")

        #Use SeedPremIDs to update premises infection status
        if ndims(SeedPremIDs) == 0 #Only a single premises status needs updating
            PremStatus[SeedPremIDs] = TimeStepVal  #Update current disease status value
            PremHasHadInfFlag[SeedPremIDs] = 1 #Update record of infection ever being present
        elseif ndims(SeedPremIDs) == 1 #Multiple premises statuses to be updated
            PremStatus[SeedPremIDs] .= TimeStepVal
            PremHasHadInfFlag[SeedPremIDs] .= 1
            println("PremStatus[SeedPremIDs]: $(PremStatus[SeedPremIDs]))")
        else
            error("SeedPremIDs should have dimension 0 or 1. Returns dimension $(ndims(SeedPremIDs)).")
        end

        #-----------------------------------------------------------------------
        # IF APPLICABLE, ASSIGN VACCINATION STAGE TO PREMISES
        # & INITIALISE SPECIFIED PROPORTION OF PREMISES AS VACCINATED IF APPROPRIATE
        #-----------------------------------------------------------------------
        if (RunControlsFn == CullIPsAndDistanceBasedCattleVaccFn!) ||
            (RunControlsFn == CullIPsAndDistanceBasedVaccFn!)

            # Disaggregate vaccination related control parameters
            VaccEff = ControlParamVals[1]::Float64
            vacc_at_simn_start_propn = ControlParamVals[2]::Float64
            vacc_response_to_risk_propn = ControlParamVals[3]::Float64

            # Error check
            if (vacc_at_simn_start_propn + vacc_response_to_risk_propn) > 1
                error("vacc_at_simn_start_propn + vacc_response_to_risk_propn exceeds 1. Invalid.")
            end

            # Get number of premises that could vaccinate at each stage
            n_vacc_at_simn_start::Int64 = round(Int64,vacc_at_simn_start_propn*PremNum)
            n_vacc_response_to_risk::Int64 = round(Int64,vacc_response_to_risk_propn*PremNum)

            # Error check, n_not_vaccinated should not be negative
            # If even split between n_vacc_at_simn_start & n_vacc_response_to_risk,
            # AND PremNum is odd, n_vacc_at_simn_start + n_vacc_response_to_risk
            # could exceed PremNum by 1.
            if n_vacc_at_simn_start + n_vacc_response_to_risk == PremNum + 1
                n_vacc_response_to_risk = n_vacc_response_to_risk -1
            end

            # Assign number in never vaccinate state
            n_not_vaccinated::Int64 = PremNum - n_vacc_at_simn_start - n_vacc_response_to_risk
            if n_not_vaccinated < 0
                error("n_not_vaccinated is negative (= $n_not_vaccinated). Invalid.")
            end
            println("n_vacc_at_simn_start: $n_vacc_at_simn_start; n_vacc_response_to_risk: $n_vacc_response_to_risk; n_not_vaccinated: $n_not_vaccinated")

            # Construct list of integers, number of entries in class X matches size of class X
            # 0 - Not ever vaccinated
            # 1 - Vaccinated prior to simulation start
            # 2 - Vaccinated when risk based measure is satisfied
            integer_list = [zeros(Int64,n_not_vaccinated);
                            ones(Int64,n_vacc_at_simn_start);
                            2*ones(Int64,n_vacc_response_to_risk)]

            # Update vector tracking the stage each premises vaccinates at
            # 0 - Not ever vaccinated
            # 1 - Vaccinated prior to simulation start
            # 2 - Vaccinated when risk based measure is satisfied
            ControlParamVals[5]::Vector{Int64} .= sample(rng,integer_list, PremNum, replace = false)
                # Done by sampling from list of integers without repetition

            # Sample the inoculation times for each premises
            # Time between vaccination being administered and becoming effective (if not infected)
            prem_time_to_inoculation_fn::Function = ControlParamVals[6]
            ControlParamVals[7] .= prem_time_to_inoculation_fn(rng,PremNum)

            # For those assigned as being vaccinated prior to simulation start
            # update vaccination related variables (if not assigned as a seed infected)
            # Only cattle are vaccinated
            for prem_idx = 1:PremNum
                if (0<=PremStatus[prem_idx]<DetectionTime) &&
                        (ControlParamVals[5][prem_idx] == 1) &&
                        (prem_idx ∉ SeedPremIDs)
                    # Premises to be vaccinated from start of simn
                    # AND premises is not a seed infected (if seed infected, assume did not have vaccine)

                    # Vaccinated premises, update PremVaccStatus
                    PremVaccStatus[prem_idx] = 1

                    # Amend prem_time_to_inoculation to be zero
                    # Corresponds to ControlParamVals[7][prem_idx]
                    ControlParamVals[7][prem_idx] = 0

                    if (RunControlsFn == CullIPsAndDistanceBasedVaccFn!)
                        # All livestock types are vaccinated

                        # Update SpeciesGroupVaccStatusByPrem
                        SpeciesGroupVaccStatusByPrem[prem_idx,:] .= 1

                        # Check premises is susceptible & not a seed infected
                        # If premises is a seed infected, vaccination is ineffective
                        if (PremStatus[prem_idx] == 0)
                            PremSuscept_SpeciesBreakdown[prem_idx,:] = PremSuscept_SpeciesBreakdown[prem_idx,:].*(1-VaccEff)
                            PremTransmiss_SpeciesBreakdown[prem_idx,:] = PremTransmiss_SpeciesBreakdown[prem_idx,:].*(1-VaccEff)
                            PremSuscept[prem_idx] = sum(PremSuscept_SpeciesBreakdown[prem_idx,:])
                            PremTransmiss[prem_idx] = sum(PremTransmiss_SpeciesBreakdown[prem_idx,:])
                        end
                    else
                        # Only cattle are vaccinated (first column of livestock type)

                        # Update SpeciesGroupVaccStatusByPrem
                        SpeciesGroupVaccStatusByPrem[prem_idx,1] = 1

                        # Check premises is susceptible & not a seed infected
                        # If premises is a seed infected, vaccination is ineffective
                        if (PremStatus[prem_idx] == 0)
                            PremSuscept_SpeciesBreakdown[prem_idx,1] = PremSuscept_SpeciesBreakdown[prem_idx,1].*(1-VaccEff)
                            PremTransmiss_SpeciesBreakdown[prem_idx,1] = PremTransmiss_SpeciesBreakdown[prem_idx,1].*(1-VaccEff)
                            PremSuscept[prem_idx] = sum(PremSuscept_SpeciesBreakdown[prem_idx,:])
                            PremTransmiss[prem_idx] = sum(PremTransmiss_SpeciesBreakdown[prem_idx,:])
                        end
                    end
                end
            end
        end

        #-----------------------------------------------------------------------
        # INITIALISE PREMISES SUSCEPTIBILITIES
        #-----------------------------------------------------------------------
        # - Drawing susceptibilities relative to seeding of RNG

        # Reinitialise the RNG
        rng = MersenneTwister(RNGseed+ItrIdx)

        # Generate a random number according to the exponential distribution with scale 1.
        prem_suscep_remaining_sellke = randexp(rng,PremNum)

        #-----------------------------------------------------------------------
        # RUN REPLICATE
        #-----------------------------------------------------------------------
        #Profile.clear() #Clear profiler
        @time RunSimulationReplicateFn(BatchID,
                                        ItrIdx,
                                        EpiParamVals,
                                        ControlParamVals,
                                        LinkedPremToControlIdxs,
                                        CoordType,
                                        PremStatus,
                                        PremHasHadInfFlag,
                                        PremVaccStatus,
                                        SpeciesGroupVaccStatusByPrem,
                                        EventIndicatorArray,
                                        EventTimeArray,
                                        PremLoc_AllVals,
                                        PremLoc_xVals,
                                        PremLoc_yVals,
                                        PremLivestockData,
                                        PremSuscept,
                                        PremTransmiss,
                                        PremSuscept_SpeciesBreakdown,
                                        PremTransmiss_SpeciesBreakdown,
                                        OutputLevels,
                                        MaxTime,
                                        TimeStepVal,
                                        KernelLookUpVec,
                                        RunControlsFn,
                                        ReplicateFilePrefix,
                                        OutputFileObjs,
                                        IterateOutbreakFn,
                                        prem_suscep_remaining_sellke,
                                        precalc_dist_and_kernel_array_flag,
                                        kernel_val_array)
        #Profile.print(format=:flat)
    end

    #Close files (for files storing data across all nReps)
    for FileIdx = 1:length(OutputFileObjs)
        close(OutputFileObjs[FileIdx])
    end

    return nothing

end

#-------------------------------------------------------------------------------
# LOAD & RUN DESIRED CONFIGURATION
#-------------------------------------------------------------------------------
"""
    RunSpatialSellkeSimn(SellkeConfigFn::Function,
                            BatchID::String,
                            nReps::Int64,
                            RNGseed::Int64,
                            precalc_dist_and_kernel_array_flag::Bool,
                            control_param_input_args::Array{Float64,1},
                            config_txt_file_name::String,
                            config_JLD2_file_name::String)

Load and run spatial simulation on selected configuration using the Sellke construction.

Inputs:
- `SellkeConfigFn::Function`: Load desired variable configuration from a premade function
- `BatchID::String`: String to be used in filenames as an identifier
- `nReps::Int64`: Number of replicates to run using selected configuration
- `RNGseed::Int64`: Input to seed RNG
- `precalc_dist_and_kernel_array_flag::Bool`: Indicator. If true, construct array of premises-to-premises kernel quantities before running main outbreak function
- `control_param_input_args::Array{Float64,1}`:  Vector: [Propn controlled prior to start time, Propn controlled if risk threshold surpassed (only control if still eligible), Risk measure]
- `config_txt_file_name::String`: Text file to keep record of variable values used in each BatchID
- `config_JLD2_file_name::String`: JLD2 file to keep record of variable values used in each BatchID

Outputs: None \n
Location: run\\_generic\\_livestock\\_disease\\_spatial\\_model\\_sellke\\_vers.jl
"""
function RunSpatialSellkeSimn(SellkeConfigFn::Function,
                                BatchID::String,
                                nReps::Int64,
                                RNGseed::Int64,
                                precalc_dist_and_kernel_array_flag::Bool,
                                control_param_input_args::Array{Float64,1},
                                config_txt_file_name::String,
                                config_JLD2_file_name::String)

    #Initialise RNG seed
    Random.seed!(RNGseed)
    rng = MersenneTwister(RNGseed)

    #Load variable configuration
    PremLocData, PremLivestockData, TimeParams, CalcPremSuscepTransmissFn, KernelFn,
            EpiParamVals,
            ControlParamVals,
            LinkedPremToControlIdxs,
            InitialInfInfo,
            OutputLevels,
            PerAnimalSuscep,
            PerAnimalTransmiss,
            SuscepExponent,
            TransmissExponent,
            IterateOutbreakFn,
            RunControlsFn,
            RunSimulationReplicateFn,
            ReplicateFilePrefix,
            OutputFileObjs = SellkeConfigFn(rng,control_param_input_args)

    println("Loaded variable configuration.")

    # Save config details to txt file
    io = open(config_txt_file_name, "w")
    writedlm(io, ["ConfigFn: $SellkeConfigFn"])
    writedlm(io, ["BatchID: $BatchID"])
    writedlm(io, ["nReps: $nReps"])
    writedlm(io, ["RNGseed: $RNGseed"])
    writedlm(io, ["precalc_dist_and_kernel_array_flag: $precalc_dist_and_kernel_array_flag"])
    writedlm(io, ["TimeParams: $TimeParams"])
    writedlm(io, ["CalcPremSuscepTransmissFn: $CalcPremSuscepTransmissFn"])
    writedlm(io, ["KernelFn: $KernelFn"])
    writedlm(io, ["EpiParamVals: $EpiParamVals"])
    writedlm(io, ["ControlParamVals: $ControlParamVals"])
    writedlm(io, ["InitialInfInfo: $InitialInfInfo"])
    writedlm(io, ["OutputLevels: $OutputLevels"])
    writedlm(io, ["PerAnimalSuscep: $PerAnimalSuscep"])
    writedlm(io, ["PerAnimalTransmiss: $PerAnimalTransmiss"])
    writedlm(io, ["SuscepExponent: $SuscepExponent"])
    writedlm(io, ["TransmissExponent: $TransmissExponent"])
    writedlm(io, ["IterateOutbreakFn: $IterateOutbreakFn"])
    writedlm(io, ["RunControlsFn: $RunControlsFn"])
    writedlm(io, ["RunSimulationReplicateFn: $RunSimulationReplicateFn"])
    writedlm(io, ["ReplicateFilePrefix: $ReplicateFilePrefix"])
    close(io)

    # Save config details to JLD 2 file
    jldsave(config_JLD2_file_name;
            SellkeConfigFn,
            BatchID,
            nReps,
            RNGseed,
            PremLocData,
            PremLivestockData,
            TimeParams,
            CalcPremSuscepTransmissFn,
            KernelFn,
            EpiParamVals,
            ControlParamVals,
            LinkedPremToControlIdxs,
            InitialInfInfo,
            OutputLevels,
            PerAnimalSuscep,
            PerAnimalTransmiss,
            SuscepExponent,
            TransmissExponent,
            IterateOutbreakFn,
            RunControlsFn,
            RunSimulationReplicateFn,
            ReplicateFilePrefix,
            OutputFileObjs)

    #Call model simulation function (that uses the Sellke construction)
    SpatialLivestockSellkeInfectionModel(BatchID,PremLocData,PremLivestockData::Array{Int64},
                                                        TimeParams,
                                                        nReps::Int64,
                                                        CalcPremSuscepTransmissFn,
                                                        KernelFn,
                                                        EpiParamVals,
                                                        ControlParamVals,
                                                        LinkedPremToControlIdxs,
                                                        InitialInfInfo,
                                                        OutputLevels,
                                                        PerAnimalSuscep::Array{Float64},
                                                        PerAnimalTransmiss::Array{Float64},
                                                        SuscepExponent::Array{Float64},
                                                        TransmissExponent::Array{Float64},
                                                        IterateOutbreakFn::Function,
                                                        RunControlsFn,
                                                        RunSimulationReplicateFn::Function,
                                                        ReplicateFilePrefix::String,
                                                        OutputFileObjs::Array{IOStream,1},
                                                        precalc_dist_and_kernel_array_flag::Bool)
    return nothing
end

#-------------------------------------------------------------------------------
# SET VARIABLES FROM ARGS
#-------------------------------------------------------------------------------
args = ARGS

# If running locally from REPL, will not have any input ARGS
# Can set values for input parameters to main run function here
# args/ARGS list
# args[1] job_ID
# args[2] job_ID_foundation_value: Set value that job_ID will be incremented upon
# args[3] BatchID_offset: Value that BatchID is offset by
# args[4] RNGseed: To be used to initialise the random number generator
# args[5] ConfigFn: Location and pathogen configuration
# args[6] nReps: Number of replicates requested
# args[7] precalc_dist_and_kernel_array_flag: Set if should precalculate distance and kernel arrays
if length(ARGS)==0
    args = [ "11", "0", "0", "99", "CumbriaSellkeFMDconfig", "10", "true"]
end

# To run from command line, example:
# julia run_generic_livestock_disease_spatial_model_sellke_vers.jl 1 3000 5000 99 CumbriaFMDconfig 10 true

# Set identifier for job
const job_ID = parse(Int64, args[1])

# Set value that job_ID will be incremented upon
const job_ID_foundation_value = parse(Int64, args[2])

# Set value that BatchID will be offset by
const BatchID_offset = parse(Int64, args[3])

# Set RNG seed
const RNGseed = parse(Int64, args[4])

# # Specify configuration function
# # Relevant files in "../GridSimnFns/SimnVarConfigs.jl"
s = Symbol(args[5])
const SellkeConfigFn = getfield(Main, s) #Make Symbol a callable function

# Set number of replicates to be run
const nReps = parse(Int64, args[6])

# Set if should precalculate distance and kernel arrays
const precalc_dist_and_kernel_array_flag = parse(Bool, args[7])

# Assign BatchID from list of possible BatchIDs
const BatchID_vec = string.(BatchID_offset .+ collect(1:50000))
const BatchID = BatchID_vec[job_ID_foundation_value+job_ID]


# Load control strategy variables based on the job_ID
# control_param_input_args - Vector: [Propn controlled prior to start time,
#                                     Propn controlled if risk threshold surpassed, (only control if still eligible)
#                                     Risk measure]
control_param_input_args_array = readdlm("parameter_combination_files/generic_model_control_param_val_array.txt") # Load the parameter combinations from file
control_param_input_args = control_param_input_args_array[job_ID_foundation_value+job_ID,:]


# File to store variable values
if (SellkeConfigFn == CumbriaSellkeFMDconfig) || (SellkeConfigFn == CumbriaSellkeFMDconfig_alternate_transmiss_params)
    config_txt_file_name = string("config_log_files/cumbria_models_sellke_simn_config_BatchID",BatchID,".txt")
    config_JLD2_file_name = string("config_JLD2_files/cumbria_models_grid_simn_config_BatchID",BatchID,".jld2")
elseif (SellkeConfigFn == DevonSellkeFMDconfig) || (SellkeConfigFn == DevonSellkeFMDconfig_alternate_transmiss_params)
    config_txt_file_name = string("config_log_files/devon_models_sellke_simn_config_BatchID",BatchID,".txt")
    config_JLD2_file_name = string("config_JLD2_files/devon_models_grid_simn_config_BatchID",BatchID,".jld2")
elseif (SellkeConfigFn == TestSellkeSimnConfig)
    config_txt_file_name = string("config_log_files/dummy_model_sellke_simn_config_BatchID",BatchID,".txt")
    config_JLD2_file_name = string("config_JLD2_files/dummy_model_grid_simn_config_BatchID",BatchID,".jld2")
else
    error("Invalid SellkeConfigFn provided.")
end

RunSpatialSellkeSimn(SellkeConfigFn,
                    BatchID,
                    nReps,
                    RNGseed,
                    precalc_dist_and_kernel_array_flag,
                    control_param_input_args,
                    config_txt_file_name,
                    config_JLD2_file_name)
