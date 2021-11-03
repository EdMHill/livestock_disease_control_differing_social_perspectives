#=
Purpose:
File to house variable configurations to input to spatial transmission model simulation
that uses the Sellke construction

Configuration list:
- Generic grid, FMD-like disease
- Cumbria, FMD-like disease
- Cumbria, FMD-like disease with lower transmissibility & longer infectious period
- Devon, FMD-like disease
- Devon, FMD-like disease with lower transmissibility & longer infectious period

Date: 3rd November 2021
=#

# ------------------------------------------------------------------------------
### Generic grid, FMD-like disease
# ------------------------------------------------------------------------------
"""
    TestSellkeSimnConfig(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

Setup configuration to be run using the Sellke construction: Test function to run an FMD-like disease on a generic landscape.

Inputs:
- `rng::AbstractRNG`: The random number generator.
- `control_param_input_args::Array{Float64,1}`:  Vector: [Propn controlled prior to start time, Propn controlled if risk threshold surpassed (only control if still eligible), Risk measure]

Outputs:
- `PremLocData::Array{Any,1}`: First entry: Co-ordinate type; Second entry: Columns with x/y or lat/long data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
- `PremLivestockData::Array{Int64}`: Number of each livestock type per premises. Row per premises, column per animal.
- `TimeParams::Array{Float64,1}`: Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over.
- `CalcPremSuscepTransmissFn::Function`: Use the specified function to calculate premises-level susceptibility and transmissibility
- `KernelFn::Function`: Use the specified kernel form. Defines risk of transmission w.r.t. distance.
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
- `ReplicateFilePrefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `OutputFileObjs::Array{IOStream,1}`: Filename identifiers written to by all replicates.

Location: SellkeSimnVarConfigs.jl
"""
function TestSellkeSimnConfig(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

    ### Load dummy livestock data ###
    dummy_data_filename = "../../../data/dummy_data/livestock_dummy_data.txt"
    CountyRawData = readdlm(dummy_data_filename,Int64)

    #Column breakdown
    ## Cols 1-4: Survey year, CPH  Easting, Northing,
    # Cols 5-10: Cattle, Sheep, Horses, Pigs, Geese, Chickens,
    # Cols 11-15 Ducks, Goats, Turkeys, Other Poultry, Deer,
    # Cols 16-20 Other livestock, Alpacas, Llamas, Bee hives

    #Retain only Cattle and Sheep (columns 5 and 6)
    cattle_and_sheep_column_idxs = [5,6]

    # Column indexes for location data
    location_column_idxs = [3,4]

    #Cast as integer type array
    PremLivestockDataAll = convert(Array{Int64,2},abs.(CountyRawData[:,cattle_and_sheep_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Retain those premises that have cattle and/or sheep present
    TotalCattleSheepOnPrem = sum(PremLivestockDataAll,dims=2)
    RetainPremFlag = TotalCattleSheepOnPrem[:].>0 #Use [:} to ensure RetainPremFlag is type BitArray{1}, not BitArray{2}!
    PremLivestockData = PremLivestockDataAll[RetainPremFlag,:]

    #Assign location info to variable
    PremLocRaw = convert(Array{Float64,2},abs.(CountyRawData[RetainPremFlag,location_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Get premises location data
    PremNum = 5000
    PremLocRaw = zeros(PremNum,2)
    xVals = rand(rng,PremNum).*100000
    yVals = rand(rng,PremNum).*100000
    PremLocRaw = [xVals yVals]

    #Get premises location data
    # PremLocData - (tuple) First entry: Co-ordinate type; Second entry: Columns with x/y or long/lat data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
    BoundingBoxVar = [0,100000.,0.,100000.]
    PremLocData = [1, #CoordType setting, 3 for LatLong
                    PremLocRaw,
                    BoundingBoxVar]

    # PremLivestockData - (integer, array)  Number of each livestock type per premises. Row per premises, column per animal.
    PremLivestockData = rand(100:1:1000,PremNum,2)

    #-------------------------------------------------------------------------------
    ### ERROR CHECKS
    #-------------------------------------------------------------------------------
    #CHECK THAT ROW COUNTS MATCH WHEN COMPARING PREMISES LOCATION DATA & LIVESTOCK POPN DATA!
    PremLocRaw = PremLocData[2]
    if size(PremLocRaw,1) != size(PremLivestockData,1)
        error("Inconsistency in number of records in location dataset($(size(PremLocDataCoords,1))) and number of records in livestock dataset, $(size(PremLivestockData,1)).")
    end

    CheckLandscapeValid(BoundingBoxVar,
                                PremLocRaw[:,1],
                                PremLocRaw[:,2])
    #-------------------------------------------------------------------------------

    # TimeParams - Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over
    TimeStepVal = 1.0; MaxTime = 10*365.;
    TimeParams = [TimeStepVal,MaxTime]

    # CalcPremSuscepTransmissFn - Use the specified function to calculate premises-level susceptibility and transmissibility
    CalcPremSuscepTransmissFn = CalcPremSuscepTransmiss_FMDlike

    # KernelFn - (function name) Use the specified kernel form. Defines risk of transmission w.r.t. distance
    KernelFn = Construct_USDOS2_kernel

    # EpiParamVals - (tuple) [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)]
    IncubationTime = 5.; DetectionTime = 9.; RemovalTime = 13.;
    EpiParamVals = [IncubationTime, DetectionTime, RemovalTime]

    # ControlParamVals - Variables related to implementing control measures
    ControlParamVals = [1.] #Set ring cull distance

    # LinkedPremToControlIdxs - For each premises, those linked premises that would also undergo control
    LinkedPremToControlIdxs = Array{Array{Int64,1}}(undef)

    # InitialInfInfo - (tuple) [SeedMethod,NumOfNodes/NodeIDs]
    #                              Seed method. 1 = random,
    #                                            2 = single specific node id,
    #                                            3 = set of specific node ids,
    #                                            4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
    #                                           5 = seed a random site and it's N nearest neighbours
    #                              Number of nodes to seed each replicate if SeedMethod = 1 or 5; or if SeedMethod = 2/3, seed this specific node every replicate.
    #                                   Node ids to seed. One id/row, number of lines must be == number of replicates.
    SeedMethod = 5
    NumOfNodes_Or_NodeIDs = 5
    InitialInfInfo = [SeedMethod,NumOfNodes_Or_NodeIDs]

    # OutputLevels - (Vector) Will output replicate results when the cumulative number of infected reaches these levels. Will also exit after the last output. Set to 0 if want to output every timestep.
    OutputLevels = 0

    # PerAnimalSuscep -  Susceptibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalSuscep = [5.7 1]

    # PerAnimalTransmiss -  Transmissibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalTransmiss = [8.2e-4 8.3e-4]

    # SuscepExponent -  Susceptibility exponents for each livestock type
    # [Cattle, Sheep]
    SuscepExponent = [0.41 0.2]

    # TransmissExponent -  Transmissibility exponents for each livestock type
    # [Cattle, Sheep]
    TransmissExponent = [0.42 0.49]

    # IterateOutbreakFn - (function) Function to perform disease transitions and control implementation per timestep
    IterateOutbreakFn = IterateOutbreak_TestSellkeFMDsimn!

    #  RunControlsFn - Enact controls as specified within the given function.
    RunControlsFn = CullIPsAndDistanceBasedVaccFn!
    PremVaccStage = zeros(Int64,PremNum)

    # Set up delay in vaccine becoming effective. Assign a value per premises
    prem_time_to_inoculation_fn = uniform_four_to_six_days_inoculation_fn
    prem_time_to_inoculation_vec = zeros(Float64,PremNum) # Placeholder vector of inoculation time for each presmies. Written to during each replicate

    # Aggregate control paramter variables
    ControlParamVals = [1.,
                        control_param_input_args[1],
                        control_param_input_args[2],
                        control_param_input_args[3],
                        PremVaccStage,
                        prem_time_to_inoculation_fn,
                        prem_time_to_inoculation_vec]
        #[Vacc efficacy,
        #    Proportion vaccinated prior to start time,
        #    Proportion vaccinated if risk threshold surpassed, (only vaccinate if still eligible)
        #   Distance threshold (in metres): Newly notified infection within this range can trigger vaccination
        #   vector to specify when a premises would undergo vaccination,
        #   Function for drawing inoculation times,
        #   vector to specify length of time for vaccination to become effective]

    # RunSimulationReplicateFn - (function) Function to perform single outbreak replicate
    RunSimulationReplicateFn = RunTestSellkeFMDsimn

    #ReplicateFilePrefix - Directory location to be used with output files storing data with individual file per replicate
    ReplicateFilePrefix = "../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs"

    #  OutputFileObjs - Filename identifiers. Used for files written to by all replicates.
    OutbreakDurationFile = open("../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs_Aggregated/OutbreakDuration_BatchID$(BatchID).txt", "a")
    PremPerDiseaseStateFile = open("../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs_Aggregated/PremPerDiseaseState_BatchID$(BatchID).txt", "a")
    CumulativeCulledFile = open("../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs_Aggregated/CumulativeCulled_BatchID$(BatchID).txt", "a")
    CumulativeVaccFile = open("../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs_Aggregated/CumulativeVacc_BatchID$(BatchID).txt", "a")
    CumulativeCasesAnimalsFile = open("../../../results/GenericLandscapeSimnOutputs/GenericLandscape_EpiOutputs_Aggregated/CumulativeCasesAnimals_BatchID$(BatchID).txt", "a")
    OutputFileObjs = [OutbreakDurationFile,PremPerDiseaseStateFile,CumulativeCulledFile,CumulativeVaccFile,CumulativeCasesAnimalsFile]

    return PremLocData::Array{Any,1},
            PremLivestockData::Array{Int64},
            TimeParams::Array{Float64,1},
            CalcPremSuscepTransmissFn::Function,
            KernelFn::Function,
            EpiParamVals::Array{Float64,1},
            ControlParamVals::Array{Any,1},
            LinkedPremToControlIdxs::Array{Array{Int64,1}},
            InitialInfInfo,
            OutputLevels,
            PerAnimalSuscep::Array{Float64,2},
            PerAnimalTransmiss::Array{Float64,2},
            SuscepExponent::Array{Float64,2},
            TransmissExponent::Array{Float64,2},
            IterateOutbreakFn::Function,
            RunControlsFn::Function,
            RunSimulationReplicateFn::Function,
            ReplicateFilePrefix::String,
            OutputFileObjs::Array{IOStream,1}
end


#-------------------------------------------------------------------------------
### Cumbria, FMD-like disease
#-------------------------------------------------------------------------------
"""
    CumbriaSellkeFMDconfig(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

Setup configuration to be run using the Sellke construction: Cumbria, FMD-like disease.

Inputs:
- `rng::AbstractRNG`: The random number generator.
- `control_param_input_args::Array{Float64,1}`:  Vector: [Propn controlled prior to start time, Propn controlled if risk threshold surpassed (only control if still eligible), Risk measure]

Outputs:
- `PremLocData::Array{Any,1}`: First entry: Co-ordinate type; Second entry: Columns with x/y or lat/long data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
- `PremLivestockData::Array{Int64}`: Number of each livestock type per premises. Row per premises, column per animal.
- `TimeParams::Array{Float64,1}`: Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over.
- `CalcPremSuscepTransmissFn::Function`: Use the specified function to calculate premises-level susceptibility and transmissibility
- `KernelFn::Function`: Use the specified kernel form. Defines risk of transmission w.r.t. distance.
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
- `ReplicateFilePrefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `OutputFileObjs::Array{IOStream,1}`: Filename identifiers written to by all replicates.

Location: SellkeSimnVarConfigs.jl
"""
function CumbriaSellkeFMDconfig(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

    ### Load Cumbria livestock data ###
    cumbria_data_filename = "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID08.txt"
    CountyRawData = readdlm(cumbria_data_filename)

    #Assign cattle herd and sheep flock size to variable
    # AND assign location info to variable
    if cumbria_data_filename == "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID08.txt"
        #Column breakdown
        ## Cols 1-4: Survey year, CPH  Easting, Northing,
        # Cols 5-10: Cattle, Sheep, Horses, Pigs, Geese, Chickens,
        # Cols 11-15 Ducks, Goats, Turkeys, Other Poultry, Deer,
        # Cols 16-20 Other livestock, Alpacas, Llamas, Bee hives

        #Retain only Cattle and Sheep (columns 5 and 6)
        cattle_and_sheep_column_idxs = [5,6]

        # Column indexes for location data
        location_column_idxs = [3,4]
    elseif cumbria_data_filename == "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_CountyID8.txt"
        #Columns 7-11 corrspond to Cattle, Pigs, Sheep, Goats, Deer.
        #Retain only Cattle and Sheep (columns 7 and 9)
        cattle_and_sheep_column_idxs = [7,9]

        # Column indexes for location data
        location_column_idxs = [4,5]
    else
        error("Unrecognised cumbria_data_filename")
    end

    #Cast as integer type array
    PremLivestockDataAll = convert(Array{Int64,2},abs.(CountyRawData[:,cattle_and_sheep_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Retain those premises that have cattle and/or sheep present
    TotalCattleSheepOnPrem = sum(PremLivestockDataAll,dims=2)
    RetainPremFlag = TotalCattleSheepOnPrem[:].>0 #Use [:} to ensure RetainPremFlag is type BitArray{1}, not BitArray{2}!
    PremLivestockData = PremLivestockDataAll[RetainPremFlag,:]

    #Assign location info to variable
    PremLocRaw = convert(Array{Float64,2},abs.(CountyRawData[RetainPremFlag,location_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Set bounding box manually!
    BoundingBoxLeft = 293400.
    BoundingBoxRight = 389900.
    BoundingBoxBottom = 460500.
    BoundingBoxTop = 595000.
    BoundingBoxVar = [BoundingBoxLeft,BoundingBoxRight,BoundingBoxBottom,BoundingBoxTop]

    # PremLocData - (tuple) First entry: Co-ordinate type; Second entry: Columns with x/y or long/lat data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
    PremLocData = [1, #"Cartesian" (distance between co-ords in metres)
                    PremLocRaw,
                    BoundingBoxVar]
    #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    ### ERROR CHECKS
    #-------------------------------------------------------------------------------
    #CHECK THAT ROW COUNTS MATCH WHEN COMPARING PREMISES LOCATION DATA & LIVESTOCK POPN DATA!
    if size(PremLocRaw,1) != size(PremLivestockData,1)
        error("Inconsistency in number of records in location dataset($(size(PremLocDataCoords,1))) and number of records in livestock dataset, $(size(PremLivestockData,1)).")
    end

    CheckLandscapeValid(BoundingBoxVar,
                                PremLocRaw[:,1],
                                PremLocRaw[:,2])
    #-------------------------------------------------------------------------------

    # TimeParams - Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over
    TimeStepVal = 1.; MaxTime = 10*365.;
    TimeParams = [TimeStepVal,MaxTime]

    # CalcPremSuscepTransmissFn - Use the specified function to calculate premises-level susceptibility and transmissibility
    CalcPremSuscepTransmissFn = CalcPremSuscepTransmiss_FMDlike

    # KernelFn - (function name) Use the specified kernel form. Defines risk of transmission w.r.t. distance
    #KernelFn = Construct_GB_FMD_kernel::Function
    KernelFn = Construct_USDOS2_kernel

    # EpiParamVals - (tuple) [Incubation time, Detection time, Removal time (culled)]
    IncubationTime = 5.; DetectionTime = 9.; RemovalTime = 13.;
    EpiParamVals = [IncubationTime, DetectionTime, RemovalTime]

    # LinkedPremToControlIdxs - For each premises, those linked premises that would also undergo control
    PremNum = size(PremLocRaw,1) #Get number of premises
    LinkedPremToControlIdxs = Array{Array{Int64,1}}(undef,PremNum) #Intialise storage tuple


    # InitialInfInfo - (tuple) [SeedMethod,NumOfNodes/NodeIDs]
    #                              Seed method. 1 = random,
    #                                            2 = single specific node id,
    #                                            3 = set of specific node ids,
    #                                            4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
    #                                           5 = seed a random site and it's N nearest neighbours
    #                              Number of nodes to seed each replicate if SeedMethod = 1 or 5; or if SeedMethod = 2/3, seed this specific node every replicate.
    #                                   Node ids to seed. One id/row, number of lines must be == number of replicates.
    SeedMethod = 5; NumOfNodes_Or_NodeIDs = 5;
    #SeedMethod = 2
    #NumOfNodes_Or_NodeIDs = 17816
    #NumOfNodes_Or_NodeIDs = [17000,17100,17200,17300,17400,17590,17789,17807,17808,17816]
    InitialInfInfo = [SeedMethod,NumOfNodes_Or_NodeIDs]

    # OutputLevels - (Vector) Will output replicate results when the cumulative number of infected reaches these levels. Will also exit after the last output. Set to 0 if want to output every timestep.
    OutputLevels = 0
    #OutputLevels = [10,100,1000]

    # PerAnimalSuscep -  Susceptibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalSuscep = [5.7 1] #From 2008 FMD paper for Cumbria

    # PerAnimalTransmiss -  Transmissibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalTransmiss = [8.2e-4 8.3e-4] #From 2008 FMD paper for Cumbria

    # SuscepExponent -  Susceptibility exponents for each livestock type
    # [Cattle, Sheep]
    SuscepExponent = [0.41 0.2] #From 2008 FMD paper for Cumbria

    # TransmissExponent -  Transmissibility exponents for each livestock type
    # [Cattle, Sheep]
    TransmissExponent = [0.42 0.49] #From 2008 FMD paper for Cumbria

    # IterateOutbreakFn - (function) Function to perform disease transitions and control implementation per timestep
    IterateOutbreakFn = IterateOutbreak_TestSellkeFMDsimn!

    #  RunControlsFn - Enact controls as specified within the given function.
    # ControlParamVals - Variables related to implementing control measures
    RunControlsFn = CullIPsAndDistanceBasedVaccFn! #CullIPsAndPremPrevalenceBasedVaccFn! # CullIPsAndCumulPremCaseBasedVaccFn!
    PremVaccStage = zeros(Int64,PremNum)

    # Set up delay in vaccine becoming effective. Assign a value per premises
    prem_time_to_inoculation_fn = uniform_four_to_six_days_inoculation_fn
    prem_time_to_inoculation_vec = zeros(Float64,PremNum) # Placeholder vector of inoculation time for each presmies. Written to during each replicate

    # Aggregate control paramter variables
    ControlParamVals = [1.,
                        control_param_input_args[1],
                        control_param_input_args[2],
                        control_param_input_args[3],
                        PremVaccStage,
                        prem_time_to_inoculation_fn,
                        prem_time_to_inoculation_vec]
        #[Vacc efficacy,
        #    Proportion vaccinated prior to start time,
        #    Proportion vaccinated if risk threshold surpassed, (only vaccinate if still eligible)
        #   Distance threshold (in metres): Newly notified infection within this range can trigger vaccination
        #   vector to specify when a premises would undergo vaccination,
        #   Function for drawing inoculation times,
        #   vector to specify length of time for vaccination to become effective]

    # RunSimulationReplicateFn - (function) Function to perform single outbreak replicate
    RunSimulationReplicateFn = RunTestSellkeFMDsimn

    #ReplicateFilePrefix - Directory location to be used with output files storing data with individual file per replicate
    ReplicateFilePrefix = "../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs"

    #  OutputFileObjs - Filename identifiers. Used for files written to by all replicates.
    OutbreakDurationFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/OutbreakDuration_BatchID$(BatchID).txt", "a")
    PremPerDiseaseStateFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/PremPerDiseaseState_BatchID$(BatchID).txt", "a")
    CumulativeCulledFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/CumulativeCulled_BatchID$(BatchID).txt", "a")
    CumulativeVaccFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/CumulativeVacc_BatchID$(BatchID).txt", "a")
    CumulativeCasesAnimalsFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/CumulativeCasesAnimals_BatchID$(BatchID).txt", "a")
    OutputFileObjs = [OutbreakDurationFile,PremPerDiseaseStateFile,CumulativeCulledFile,CumulativeVaccFile,CumulativeCasesAnimalsFile]

    return PremLocData::Array{Any,1},
            PremLivestockData::Array{Int64},
            TimeParams::Array{Float64,1},
            CalcPremSuscepTransmissFn::Function,
            KernelFn::Function,
            EpiParamVals::Array{Float64,1},
            ControlParamVals::Array{Any,1},
            LinkedPremToControlIdxs::Array{Array{Int64,1}},
            InitialInfInfo,
            OutputLevels,
            PerAnimalSuscep::Array{Float64,2},
            PerAnimalTransmiss::Array{Float64,2},
            SuscepExponent::Array{Float64,2},
            TransmissExponent::Array{Float64,2},
            IterateOutbreakFn::Function,
            RunControlsFn::Function,
            RunSimulationReplicateFn::Function,
            ReplicateFilePrefix::String,
            OutputFileObjs::Array{IOStream,1}
end

#-------------------------------------------------------------------------------
### Cumbria, FMD-like disease with shorter transmissibility & longer infectious period
#-------------------------------------------------------------------------------
"""
    CumbriaSellkeFMDconfig_alternate_transmiss_params(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

Setup configuration to be run using the Sellke construction: Cumbria, FMD-like disease with shorter transmissibility & longer infectious period.

Inputs:
- `rng::AbstractRNG`: The random number generator.
- `control_param_input_args::Array{Float64,1}`:  Vector: [Propn controlled prior to start time, Propn controlled if risk threshold surpassed (only control if still eligible), Risk measure]

Outputs:
- `PremLocData::Array{Any,1}`: First entry: Co-ordinate type; Second entry: Columns with x/y or lat/long data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
- `PremLivestockData::Array{Int64}`: Number of each livestock type per premises. Row per premises, column per animal.
- `TimeParams::Array{Float64,1}`: Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over.
- `CalcPremSuscepTransmissFn::Function`: Use the specified function to calculate premises-level susceptibility and transmissibility
- `KernelFn::Function`: Use the specified kernel form. Defines risk of transmission w.r.t. distance.
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
- `ReplicateFilePrefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `OutputFileObjs::Array{IOStream,1}`: Filename identifiers written to by all replicates.

Location: SellkeSimnVarConfigs.jl
"""
function CumbriaSellkeFMDconfig_alternate_transmiss_params(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

    ### Load Cumbria livestock data ###
    cumbria_data_filename = "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID08.txt"
    CountyRawData = readdlm(cumbria_data_filename)

    #Assign cattle herd and sheep flock size to variable
    # AND assign location info to variable
    if cumbria_data_filename == "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID08.txt"
        #Column breakdown
        ## Cols 1-4: Survey year, CPH  Easting, Northing,
        # Cols 5-10: Cattle, Sheep, Horses, Pigs, Geese, Chickens,
        # Cols 11-15 Ducks, Goats, Turkeys, Other Poultry, Deer,
        # Cols 16-20 Other livestock, Alpacas, Llamas, Bee hives

        #Retain only Cattle and Sheep (columns 5 and 6)
        cattle_and_sheep_column_idxs = [5,6]

        # Column indexes for location data
        location_column_idxs = [3,4]
    elseif cumbria_data_filename == "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_CountyID8.txt"
        #Columns 7-11 corrspond to Cattle, Pigs, Sheep, Goats, Deer.
        #Retain only Cattle and Sheep (columns 7 and 9)
        cattle_and_sheep_column_idxs = [7,9]

        # Column indexes for location data
        location_column_idxs = [4,5]
    else
        error("Unrecognised cumbria_data_filename")
    end

    #Cast as integer type array
    PremLivestockDataAll = convert(Array{Int64,2},abs.(CountyRawData[:,cattle_and_sheep_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Retain those premises that have cattle and/or sheep present
    TotalCattleSheepOnPrem = sum(PremLivestockDataAll,dims=2)
    RetainPremFlag = TotalCattleSheepOnPrem[:].>0 #Use [:} to ensure RetainPremFlag is type BitArray{1}, not BitArray{2}!
    PremLivestockData = PremLivestockDataAll[RetainPremFlag,:]

    #Assign location info to variable
    PremLocRaw = convert(Array{Float64,2},abs.(CountyRawData[RetainPremFlag,location_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Set bounding box manually!
    BoundingBoxLeft = 293400.
    BoundingBoxRight = 389900.
    BoundingBoxBottom = 460500.
    BoundingBoxTop = 595000.
    BoundingBoxVar = [BoundingBoxLeft,BoundingBoxRight,BoundingBoxBottom,BoundingBoxTop]

    # PremLocData - (tuple) First entry: Co-ordinate type; Second entry: Columns with x/y or long/lat data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
    PremLocData = [1, #"Cartesian" (distance between co-ords in metres)
                    PremLocRaw,
                    BoundingBoxVar]
    #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    ### ERROR CHECKS
    #-------------------------------------------------------------------------------
    #CHECK THAT ROW COUNTS MATCH WHEN COMPARING PREMISES LOCATION DATA & LIVESTOCK POPN DATA!
    if size(PremLocRaw,1) != size(PremLivestockData,1)
        error("Inconsistency in number of records in location dataset($(size(PremLocDataCoords,1))) and number of records in livestock dataset, $(size(PremLivestockData,1)).")
    end

    CheckLandscapeValid(BoundingBoxVar,
                                PremLocRaw[:,1],
                                PremLocRaw[:,2])
    #-------------------------------------------------------------------------------

    # TimeParams - Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over
    TimeStepVal = 1.; MaxTime = 10*365.;
    TimeParams = [TimeStepVal,MaxTime]

    # CalcPremSuscepTransmissFn - Use the specified function to calculate premises-level susceptibility and transmissibility
    CalcPremSuscepTransmissFn = CalcPremSuscepTransmiss_FMDlike

    # KernelFn - (function name) Use the specified kernel form. Defines risk of transmission w.r.t. distance
    #KernelFn = Construct_GB_FMD_kernel::Function
    KernelFn = Construct_USDOS2_kernel

    # EpiParamVals - (tuple) [Incubation time, Detection time, Removal time (culled)]
    IncubationTime = 5.; DetectionTime = 9.; RemovalTime = 37.;
    EpiParamVals = [IncubationTime, DetectionTime, RemovalTime]

    # LinkedPremToControlIdxs - For each premises, those linked premises that would also undergo control
    PremNum = size(PremLocRaw,1) #Get number of premises
    LinkedPremToControlIdxs = Array{Array{Int64,1}}(undef,PremNum) #Intialise storage tuple


    # InitialInfInfo - (tuple) [SeedMethod,NumOfNodes/NodeIDs]
    #                              Seed method. 1 = random,
    #                                            2 = single specific node id,
    #                                            3 = set of specific node ids,
    #                                            4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
    #                                           5 = seed a random site and it's N nearest neighbours
    #                              Number of nodes to seed each replicate if SeedMethod = 1 or 5; or if SeedMethod = 2/3, seed this specific node every replicate.
    #                                   Node ids to seed. One id/row, number of lines must be == number of replicates.
    SeedMethod = 5; NumOfNodes_Or_NodeIDs = 5;
    InitialInfInfo = [SeedMethod,NumOfNodes_Or_NodeIDs]

    # OutputLevels - (Vector) Will output replicate results when the cumulative number of infected reaches these levels. Will also exit after the last output. Set to 0 if want to output every timestep.
    OutputLevels = 0

    # PerAnimalSuscep -  Susceptibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalSuscep = [5.7 1] #From 2008 FMD paper for Cumbria

    # PerAnimalTransmiss -  Transmissibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalTransmiss = [8.2e-4 8.3e-4].*0.25 #From 2008 FMD paper for Cumbria

    # SuscepExponent -  Susceptibility exponents for each livestock type
    # [Cattle, Sheep]
    SuscepExponent = [0.41 0.2] #From 2008 FMD paper for Cumbria

    # TransmissExponent -  Transmissibility exponents for each livestock type
    # [Cattle, Sheep]
    TransmissExponent = [0.42 0.49] #From 2008 FMD paper for Cumbria

    # IterateOutbreakFn - (function) Function to perform disease transitions and control implementation per timestep
    IterateOutbreakFn = IterateOutbreak_TestSellkeFMDsimn!

    # RunControlsFn - Enact controls as specified within the given function.
    # ControlParamVals - Variables related to implementing control measures
    RunControlsFn = CullIPsAndDistanceBasedVaccFn! #CullIPsAndPremPrevalenceBasedVaccFn! # CullIPsAndCumulPremCaseBasedVaccFn!
    PremVaccStage = zeros(Int64,PremNum)

    # Set up delay in vaccine becoming effective. Assign a value per premises
    prem_time_to_inoculation_fn = uniform_four_to_six_days_inoculation_fn
    prem_time_to_inoculation_vec = zeros(Float64,PremNum) # Placeholder vector of inoculation time for each presmies. Written to during each replicate

    # Aggregate control paramter variables
    ControlParamVals = [1.,
                        control_param_input_args[1],
                        control_param_input_args[2],
                        control_param_input_args[3],
                        PremVaccStage,
                        prem_time_to_inoculation_fn,
                        prem_time_to_inoculation_vec]
        #[Vacc efficacy,
        #    Proportion vaccinated prior to start time,
        #    Proportion vaccinated if risk threshold surpassed, (only vaccinate if still eligible)
        #   Distance threshold (in metres): Newly notified infection within this range can trigger vaccination
        #   vector to specify when a premises would undergo vaccination,
        #   Function for drawing inoculation times,
        #   vector to specify length of time for vaccination to become effective]

    # RunSimulationReplicateFn - (function) Function to perform single outbreak replicate
    RunSimulationReplicateFn = RunTestSellkeFMDsimn

    #ReplicateFilePrefix - Directory location to be used with output files storing data with individual file per replicate
    ReplicateFilePrefix = "../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs"

    #  OutputFileObjs - Filename identifiers. Used for files written to by all replicates.
    OutbreakDurationFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/OutbreakDuration_BatchID$(BatchID).txt", "a")
    PremPerDiseaseStateFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/PremPerDiseaseState_BatchID$(BatchID).txt", "a")
    CumulativeCulledFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/CumulativeCulled_BatchID$(BatchID).txt", "a")
    CumulativeVaccFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/CumulativeVacc_BatchID$(BatchID).txt", "a")
    CumulativeCasesAnimalsFile = open("../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/CumulativeCasesAnimals_BatchID$(BatchID).txt", "a")
    OutputFileObjs = [OutbreakDurationFile,PremPerDiseaseStateFile,CumulativeCulledFile,CumulativeVaccFile,CumulativeCasesAnimalsFile]

    return PremLocData::Array{Any,1},
            PremLivestockData::Array{Int64},
            TimeParams::Array{Float64,1},
            CalcPremSuscepTransmissFn::Function,
            KernelFn::Function,
            EpiParamVals::Array{Float64,1},
            ControlParamVals::Array{Any,1},
            LinkedPremToControlIdxs::Array{Array{Int64,1}},
            InitialInfInfo,
            OutputLevels,
            PerAnimalSuscep::Array{Float64,2},
            PerAnimalTransmiss::Array{Float64,2},
            SuscepExponent::Array{Float64,2},
            TransmissExponent::Array{Float64,2},
            IterateOutbreakFn::Function,
            RunControlsFn::Function,
            RunSimulationReplicateFn::Function,
            ReplicateFilePrefix::String,
            OutputFileObjs::Array{IOStream,1}
end


#-------------------------------------------------------------------------------
### Devon, FMD-like disease
#-------------------------------------------------------------------------------
function DevonSellkeFMDconfig(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

    ### Load Devon livestock data ###
    devon_data_filename = "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID10.txt"
    CountyRawData = readdlm(devon_data_filename)

    #Assign cattle herd and sheep flock size to variable
    # AND assign location info to variable
    if devon_data_filename == "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID10.txt"
        #Column breakdown
        ## Cols 1-4: Survey year, CPH  Easting, Northing,
        # Cols 5-10: Cattle, Sheep, Horses, Pigs, Geese, Chickens,
        # Cols 11-15 Ducks, Goats, Turkeys, Other Poultry, Deer,
        # Cols 16-20 Other livestock, Alpacas, Llamas, Bee hives

        #Retain only Cattle and Sheep (columns 5 and 6)
        cattle_and_sheep_column_idxs = [5,6]

        # Column indexes for location data
        location_column_idxs = [3,4]
    elseif devon_data_filename == "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_CountyID10.txt"
        #Columns 7-11 corrspond to Cattle, Pigs, Sheep, Goats, Deer.
        #Retain only Cattle and Sheep (columns 7 and 9)
        cattle_and_sheep_column_idxs = [7,9]

        # Column indexes for location data
        location_column_idxs = [4,5]
    else
        error("Unrecognised devon_data_filename")
    end

    #Cast as integer type array
    PremLivestockDataAll = convert(Array{Int64,2},abs.(CountyRawData[:,cattle_and_sheep_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Retain those premises that have cattle and/or sheep present
    TotalCattleSheepOnPrem = sum(PremLivestockDataAll,dims=2)
    RetainPremFlag = TotalCattleSheepOnPrem[:].>0 #Use [:} to ensure RetainPremFlag is type BitArray{1}, not BitArray{2}!
    PremLivestockData = PremLivestockDataAll[RetainPremFlag,:]

    #Assign location info to variable
    PremLocRaw = convert(Array{Float64,2},abs.(CountyRawData[RetainPremFlag,location_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Set bounding box manually!
    BoundingBoxLeft = 212000.
    BoundingBoxRight = 340000.
    BoundingBoxBottom = 34960.
    BoundingBoxTop = 151210.
    BoundingBoxVar = [BoundingBoxLeft,BoundingBoxRight,BoundingBoxBottom,BoundingBoxTop]
    #-------------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    ### FIND PREMISES OUTSIDE BOUNDING BOX AND REMOVE
    #---------------------------------------------------------------------------
    valid_prem_loc_flag = (PremLocRaw[:,1] .>= BoundingBoxLeft) .&
                            (PremLocRaw[:,1] .<= BoundingBoxRight) .&
                            (PremLocRaw[:,2] .>= BoundingBoxBottom) .&
                            (PremLocRaw[:,2] .<= BoundingBoxTop)

    PremLocRaw = PremLocRaw[valid_prem_loc_flag,:]
    PremLivestockData = PremLivestockData[valid_prem_loc_flag,:]

    # PremLocData - (tuple) First entry: Co-ordinate type; Second entry: Columns with x/y or long/lat data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
    PremLocData = [1, #"Cartesian" (distance between co-ords in metres)
                    PremLocRaw,
                    BoundingBoxVar]
    #---------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    ### ERROR CHECKS
    #-------------------------------------------------------------------------------
    #CHECK THAT ROW COUNTS MATCH WHEN COMPARING PREMISES LOCATION DATA & LIVESTOCK POPN DATA!
    if size(PremLocRaw,1) != size(PremLivestockData,1)
        error("Inconsistency in number of records in location dataset($(size(PremLocDataCoords,1))) and number of records in livestock dataset, $(size(PremLivestockData,1)).")
    end

    CheckLandscapeValid(BoundingBoxVar,
                                PremLocRaw[:,1],
                                PremLocRaw[:,2])
    #-------------------------------------------------------------------------------

    # TimeParams - Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over
    TimeStepVal = 1.; MaxTime = 10*365.;
    TimeParams = [TimeStepVal,MaxTime]

    # CalcPremSuscepTransmissFn - Use the specified function to calculate premises-level susceptibility and transmissibility
    CalcPremSuscepTransmissFn = CalcPremSuscepTransmiss_FMDlike

    # KernelFn - (function name) Use the specified kernel form. Defines risk of transmission w.r.t. distance
    #KernelFn = Construct_GB_FMD_kernel::Function
    KernelFn = Construct_USDOS2_kernel

    # EpiParamVals - (tuple) [Incubation time, Detection time, Removal time (culled)]
    IncubationTime = 5.; DetectionTime = 9.; RemovalTime = 13.;
    EpiParamVals = [IncubationTime, DetectionTime, RemovalTime]

    # LinkedPremToControlIdxs - For each premises, those linked premises that would also undergo control

    PremNum = size(PremLocRaw,1) #Get number of premises
    LinkedPremToControlIdxs = Array{Array{Int64,1}}(undef,PremNum) #Intialise storage tuple


    # InitialInfInfo - (tuple) [SeedMethod,NumOfNodes/NodeIDs]
    #                              Seed method. 1 = random,
    #                                            2 = single specific node id,
    #                                            3 = set of specific node ids,
    #                                            4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
    #                                           5 = seed a random site and it's N nearest neighbours
    #                              Number of nodes to seed each replicate if SeedMethod = 1 or 5; or if SeedMethod = 2/3, seed this specific node every replicate.
    #                                   Node ids to seed. One id/row, number of lines must be == number of replicates.
    SeedMethod = 5; NumOfNodes_Or_NodeIDs = 5;
    InitialInfInfo = [SeedMethod,NumOfNodes_Or_NodeIDs]

    # OutputLevels - (Vector) Will output replicate results when the cumulative number of infected reaches these levels. Will also exit after the last output. Set to 0 if want to output every timestep.
    OutputLevels = 0

    # PerAnimalSuscep -  Susceptibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalSuscep = [5.7 1] #From 2008 FMD paper for Cumbria

    # PerAnimalTransmiss -  Transmissibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalTransmiss = [8.2e-4 8.3e-4] #From 2008 FMD paper for Cumbria

    # SuscepExponent -  Susceptibility exponents for each livestock type
    # [Cattle, Sheep]
    SuscepExponent = [0.41 0.2] #From 2008 FMD paper for Cumbria

    # TransmissExponent -  Transmissibility exponents for each livestock type
    # [Cattle, Sheep]
    TransmissExponent = [0.42 0.49] #From 2008 FMD paper for Cumbria


    # IterateOutbreakFn - (function) Function to perform disease transitions and control implementation per timestep
    IterateOutbreakFn = IterateOutbreak_TestSellkeFMDsimn!

    # RunControlsFn - Enact controls as specified within the given function.
    # ControlParamVals - Variables related to implementing control measures
    RunControlsFn = CullIPsAndDistanceBasedVaccFn! #CullIPsAndPremPrevalenceBasedVaccFn! # CullIPsAndCumulPremCaseBasedVaccFn!
    PremVaccStage = zeros(Int64,PremNum)

    # Set up delay in vaccine becoming effective. Assign a value per premises
    prem_time_to_inoculation_fn = uniform_four_to_six_days_inoculation_fn
    prem_time_to_inoculation_vec = zeros(Float64,PremNum) # Placeholder vector of inoculation time for each presmies. Written to during each replicate

    # Aggregate control paramter variables
    ControlParamVals = [1.,
                        control_param_input_args[1],
                        control_param_input_args[2],
                        control_param_input_args[3],
                        PremVaccStage,
                        prem_time_to_inoculation_fn,
                        prem_time_to_inoculation_vec]
        #[Vacc efficacy,
        #    Proportion vaccinated prior to start time,
        #    Proportion vaccinated if risk threshold surpassed, (only vaccinate if still eligible)
        #   Distance threshold (in metres): Newly notified infection within this range can trigger vaccination
        #   vector to specify when a premises would undergo vaccination,
        #   Function for drawing inoculation times,
        #   vector to specify length of time for vaccination to become effective]

    # RunSimulationReplicateFn - (function) Function to perform single outbreak replicate
    RunSimulationReplicateFn = RunTestSellkeFMDsimn

    #ReplicateFilePrefix - Directory location to be used with output files storing data with individual file per replicate
    ReplicateFilePrefix = "../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs"

    #  OutputFileObjs - Filename identifiers. Used for files written to by all replicates.
    OutbreakDurationFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/OutbreakDuration_BatchID$(BatchID).txt", "a")
    PremPerDiseaseStateFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/PremPerDiseaseState_BatchID$(BatchID).txt", "a")
    CumulativeCulledFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/CumulativeCulled_BatchID$(BatchID).txt", "a")
    CumulativeVaccFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/CumulativeVacc_BatchID$(BatchID).txt", "a")
    CumulativeCasesAnimalsFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/CumulativeCasesAnimals_BatchID$(BatchID).txt", "a")
    OutputFileObjs = [OutbreakDurationFile,PremPerDiseaseStateFile,CumulativeCulledFile,CumulativeVaccFile,CumulativeCasesAnimalsFile]

    return PremLocData::Array{Any,1},
            PremLivestockData::Array{Int64},
            TimeParams::Array{Float64,1},
            CalcPremSuscepTransmissFn::Function,
            KernelFn::Function,
            EpiParamVals::Array{Float64,1},
            ControlParamVals::Array{Any,1},
            LinkedPremToControlIdxs::Array{Array{Int64,1}},
            InitialInfInfo,
            OutputLevels,
            PerAnimalSuscep::Array{Float64,2},
            PerAnimalTransmiss::Array{Float64,2},
            SuscepExponent::Array{Float64,2},
            TransmissExponent::Array{Float64,2},
            IterateOutbreakFn::Function,
            RunControlsFn::Function,
            RunSimulationReplicateFn::Function,
            ReplicateFilePrefix::String,
            OutputFileObjs::Array{IOStream,1}
end

#-------------------------------------------------------------------------------
### Devon, FMD-like disease with shorter transmissibility & longer infectious period
#-------------------------------------------------------------------------------
"""
    DevonSellkeFMDconfig_alternate_transmiss_params(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

Setup configuration to be run using the Sellke construction: Devon, FMD-like disease with shorter transmissibility & longer infectious period.

Inputs:
- `rng::AbstractRNG`: The random number generator.
- `control_param_input_args::Array{Float64,1}`:  Vector: [Propn controlled prior to start time, Propn controlled if risk threshold surpassed (only control if still eligible), Risk measure]

Outputs:
- `PremLocData::Array{Any,1}`: First entry: Co-ordinate type; Second entry: Columns with x/y or lat/long data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
- `PremLivestockData::Array{Int64}`: Number of each livestock type per premises. Row per premises, column per animal.
- `TimeParams::Array{Float64,1}`: Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over.
- `CalcPremSuscepTransmissFn::Function`: Use the specified function to calculate premises-level susceptibility and transmissibility
- `KernelFn::Function`: Use the specified kernel form. Defines risk of transmission w.r.t. distance.
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
- `ReplicateFilePrefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `OutputFileObjs::Array{IOStream,1}`: Filename identifiers written to by all replicates.

Location: SellkeSimnVarConfigs.jl
"""
function DevonSellkeFMDconfig_alternate_transmiss_params(rng::AbstractRNG,control_param_input_args::Array{Float64,1})

    ### Load Devon livestock data ###
    devon_data_filename = "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID10.txt"
    CountyRawData = readdlm(devon_data_filename)

    #Assign cattle herd and sheep flock size to variable
    # AND assign location info to variable
    if devon_data_filename == "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_2020_CountyID10.txt"
        #Column breakdown
        ## Cols 1-4: Survey year, CPH  Easting, Northing,
        # Cols 5-10: Cattle, Sheep, Horses, Pigs, Geese, Chickens,
        # Cols 11-15 Ducks, Goats, Turkeys, Other Poultry, Deer,
        # Cols 16-20 Other livestock, Alpacas, Llamas, Bee hives

        #Retain only Cattle and Sheep (columns 5 and 6)
        cattle_and_sheep_column_idxs = [5,6]

        # Column indexes for location data
        location_column_idxs = [3,4]
    elseif devon_data_filename == "../../../data/ProcessedData/GBLivestockByCounty/GB_Farm_CountyID10.txt"
        #Columns 7-11 corrspond to Cattle, Pigs, Sheep, Goats, Deer.
        #Retain only Cattle and Sheep (columns 7 and 9)
        cattle_and_sheep_column_idxs = [7,9]

        # Column indexes for location data
        location_column_idxs = [4,5]
    else
        error("Unrecognised devon_data_filename")
    end

    #Cast as integer type array
    PremLivestockDataAll = convert(Array{Int64,2},abs.(CountyRawData[:,cattle_and_sheep_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Retain those premises that have cattle and/or sheep present
    TotalCattleSheepOnPrem = sum(PremLivestockDataAll,dims=2)
    RetainPremFlag = TotalCattleSheepOnPrem[:].>0 #Use [:} to ensure RetainPremFlag is type BitArray{1}, not BitArray{2}!
    PremLivestockData = PremLivestockDataAll[RetainPremFlag,:]

    #Assign location info to variable
    PremLocRaw = convert(Array{Float64,2},abs.(CountyRawData[RetainPremFlag,location_column_idxs])) #Negative values in data designate user-input value. Take absolute value

    #Set bounding box manually!
    BoundingBoxLeft = 212000.
    BoundingBoxRight = 340000.
    BoundingBoxBottom = 34960.
    BoundingBoxTop = 151210.
    BoundingBoxVar = [BoundingBoxLeft,BoundingBoxRight,BoundingBoxBottom,BoundingBoxTop]
    #-------------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    ### FIND PREMISES OUTSIDE BOUNDING BOX AND REMOVE
    #---------------------------------------------------------------------------
    valid_prem_loc_flag = (PremLocRaw[:,1] .>= BoundingBoxLeft) .&
                            (PremLocRaw[:,1] .<= BoundingBoxRight) .&
                            (PremLocRaw[:,2] .>= BoundingBoxBottom) .&
                            (PremLocRaw[:,2] .<= BoundingBoxTop)

    PremLocRaw = PremLocRaw[valid_prem_loc_flag,:]
    PremLivestockData = PremLivestockData[valid_prem_loc_flag,:]

    # PremLocData - (tuple) First entry: Co-ordinate type; Second entry: Columns with x/y or long/lat data; Third entry; Vector with bounding box information [Min_x,Max_x,Min_y,Max_y]
    PremLocData = [1, #"Cartesian" (distance between co-ords in metres)
                    PremLocRaw,
                    BoundingBoxVar]
    #---------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    ### ERROR CHECKS
    #-------------------------------------------------------------------------------
    #CHECK THAT ROW COUNTS MATCH WHEN COMPARING PREMISES LOCATION DATA & LIVESTOCK POPN DATA!
    if size(PremLocRaw,1) != size(PremLivestockData,1)
        error("Inconsistency in number of records in location dataset($(size(PremLocDataCoords,1))) and number of records in livestock dataset, $(size(PremLivestockData,1)).")
    end

    CheckLandscapeValid(BoundingBoxVar,
                                PremLocRaw[:,1],
                                PremLocRaw[:,2])
    #-------------------------------------------------------------------------------

    # TimeParams - Entries [TimeStepVal, MaxTime]. Timestep per iteration & timeframe simulation is run over
    TimeStepVal = 1.; MaxTime = 10*365.;
    TimeParams = [TimeStepVal,MaxTime]

    # CalcPremSuscepTransmissFn - Use the specified function to calculate premises-level susceptibility and transmissibility
    CalcPremSuscepTransmissFn = CalcPremSuscepTransmiss_FMDlike

    # KernelFn - (function name) Use the specified kernel form. Defines risk of transmission w.r.t. distance
    #KernelFn = Construct_GB_FMD_kernel::Function
    KernelFn = Construct_USDOS2_kernel

    # EpiParamVals - (tuple) [Incubation time, Detection time, Removal time (culled)]
    IncubationTime = 5.; DetectionTime = 9.; RemovalTime = 37.;
    EpiParamVals = [IncubationTime, DetectionTime, RemovalTime]

    # LinkedPremToControlIdxs - For each premises, those linked premises that would also undergo control

    PremNum = size(PremLocRaw,1) #Get number of premises
    LinkedPremToControlIdxs = Array{Array{Int64,1}}(undef,PremNum) #Intialise storage tuple


    # InitialInfInfo - (tuple) [SeedMethod,NumOfNodes/NodeIDs]
    #                              Seed method. 1 = random,
    #                                            2 = single specific node id,
    #                                            3 = set of specific node ids,
    #                                            4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
    #                                           5 = seed a random site and it's N nearest neighbours
    #                              Number of nodes to seed each replicate if SeedMethod = 1 or 5; or if SeedMethod = 2/3, seed this specific node every replicate.
    #                                   Node ids to seed. One id/row, number of lines must be == number of replicates.
    SeedMethod = 5; NumOfNodes_Or_NodeIDs = 5;
    InitialInfInfo = [SeedMethod,NumOfNodes_Or_NodeIDs]

    # OutputLevels - (Vector) Will output replicate results when the cumulative number of infected reaches these levels. Will also exit after the last output. Set to 0 if want to output every timestep.
    OutputLevels = 0

    # PerAnimalSuscep -  Susceptibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalSuscep = [5.7 1] #From 2008 FMD paper for Cumbria

    # PerAnimalTransmiss -  Transmissibility scale parameters for each livestock type
    # [Cattle, Sheep]
    PerAnimalTransmiss = [8.2e-4 8.3e-4].*0.25 #From 2008 FMD paper for Cumbria

    # SuscepExponent -  Susceptibility exponents for each livestock type
    # [Cattle, Sheep]
    SuscepExponent = [0.41 0.2] #From 2008 FMD paper for Cumbria

    # TransmissExponent -  Transmissibility exponents for each livestock type
    # [Cattle, Sheep]
    TransmissExponent = [0.42 0.49] #From 2008 FMD paper for Cumbria


    # IterateOutbreakFn - (function) Function to perform disease transitions and control implementation per timestep
    IterateOutbreakFn = IterateOutbreak_TestSellkeFMDsimn!

    # RunControlsFn - Enact controls as specified within the given function.
    # ControlParamVals - Variables related to implementing control measures
    RunControlsFn = CullIPsAndDistanceBasedVaccFn! #CullIPsAndPremPrevalenceBasedVaccFn! # CullIPsAndCumulPremCaseBasedVaccFn!
    PremVaccStage = zeros(Int64,PremNum)

    # Set up delay in vaccine becoming effective. Assign a value per premises
    prem_time_to_inoculation_fn = uniform_four_to_six_days_inoculation_fn
    prem_time_to_inoculation_vec = zeros(Float64,PremNum) # Placeholder vector of inoculation time for each presmies. Written to during each replicate

    # Aggregate control paramter variables
    ControlParamVals = [1.,
                        control_param_input_args[1],
                        control_param_input_args[2],
                        control_param_input_args[3],
                        PremVaccStage,
                        prem_time_to_inoculation_fn,
                        prem_time_to_inoculation_vec]
        #[Vacc efficacy,
        #    Proportion vaccinated prior to start time,
        #    Proportion vaccinated if risk threshold surpassed, (only vaccinate if still eligible)
        #   Distance threshold (in metres): Newly notified infection within this range can trigger vaccination
        #   vector to specify when a premises would undergo vaccination,
        #   Function for drawing inoculation times,
        #   vector to specify length of time for vaccination to become effective]

    # RunSimulationReplicateFn - (function) Function to perform single outbreak replicate
    RunSimulationReplicateFn = RunTestSellkeFMDsimn

    #ReplicateFilePrefix - Directory location to be used with output files storing data with individual file per replicate
    ReplicateFilePrefix = "../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs"

    #  OutputFileObjs - Filename identifiers. Used for files written to by all replicates.
    OutbreakDurationFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/OutbreakDuration_BatchID$(BatchID).txt", "a")
    PremPerDiseaseStateFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/PremPerDiseaseState_BatchID$(BatchID).txt", "a")
    CumulativeCulledFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/CumulativeCulled_BatchID$(BatchID).txt", "a")
    CumulativeVaccFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/CumulativeVacc_BatchID$(BatchID).txt", "a")
    CumulativeCasesAnimalsFile = open("../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/CumulativeCasesAnimals_BatchID$(BatchID).txt", "a")
    OutputFileObjs = [OutbreakDurationFile,PremPerDiseaseStateFile,CumulativeCulledFile,CumulativeVaccFile,CumulativeCasesAnimalsFile]

    return PremLocData::Array{Any,1},
            PremLivestockData::Array{Int64},
            TimeParams::Array{Float64,1},
            CalcPremSuscepTransmissFn::Function,
            KernelFn::Function,
            EpiParamVals::Array{Float64,1},
            ControlParamVals::Array{Any,1},
            LinkedPremToControlIdxs::Array{Array{Int64,1}},
            InitialInfInfo,
            OutputLevels,
            PerAnimalSuscep::Array{Float64,2},
            PerAnimalTransmiss::Array{Float64,2},
            SuscepExponent::Array{Float64,2},
            TransmissExponent::Array{Float64,2},
            IterateOutbreakFn::Function,
            RunControlsFn::Function,
            RunSimulationReplicateFn::Function,
            ReplicateFilePrefix::String,
            OutputFileObjs::Array{IOStream,1}
end
