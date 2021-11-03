#=
Purpose:
Script to find optimal risk threshold for a given relative vaccination cost
Can be done for multiple threshold values to produce an optimal threshold profile

Overview:
  - Construct vaccine to infection cost ratio values to be tested
  - Import simulation data for given threshold value (that can consist of multiple simn replicates)
  - For each ratio (and single threshold), calculate population cost (per simulation replicate)
  - For each ratio (and single threshold), calculate individual farmer perspective cost (per simulation replicate)
  - For each ratio (and single threshold), take desired cost measures across simulation replicates (e.g. mean, median, percentile,...)
  - For each ratio and summary statistic measure, across all thresholds considered, find optimal threshold to minimise cost

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
using Statistics
using MAT

#-------------------------------------------------------------------------------
### SET VARIABLES FROM ARGS
#-------------------------------------------------------------------------------
args = ARGS

# If running locally from REPL, will not have any input ARGS
# Can set values for input parameters to main run function here
# args/ARGS list
# args[1] job_ID
# args[2] ConfigFn: Location and pathogen configuration
# args[3] scen_offset: Value that scen_itr is offset by
if length(ARGS)==0
    args = ["1","CumbriaSellkeFMDconfig_alternate_transmiss_params","0"]
end

# To run from command line, example:
# julia ComputeOptimalRiskThreshold_VaccBehav_Script.jl 1 CumbriaFMDconfig 5000 300

# Set identifier for job
const job_ID = parse(Int64, args[1])

# # Relevant files in "../GridSimnFns/SimnVarConfigs.jl"
const ConfigFn = args[2] #Make Symbol a callable function

# Set value that BatchID will be offset by
const scen_offset = parse(Int64, args[3])

#-------------------------------------------------------------------------------
### CONSTRUCT VACCINE TO INFECTION COST RATIO VALUES TO BE TESTED
#-------------------------------------------------------------------------------
const VaccToInfCostRatio = 0:0.01:1 #We set C_I = 1
const VaccToInfCostRatio_NumTested = length(VaccToInfCostRatio) #Assign amount of ratios tested to variable


#-------------------------------------------------------------------------------
### SPECIFY INPUT PREFIX & BatchID_offset TO READ INPUT DATA FROM
#-------------------------------------------------------------------------------
if (ConfigFn == "CumbriaSellkeFMDconfig") ||
        (ConfigFn == "CumbriaSellkeFMDconfig_alternate_transmiss_params")
    InputFilesPrefix = "../../../results/GB_county_model_simn_outputs/Cumbria_EpiOutputs_Aggregated/"

elseif (ConfigFn == "DevonSellkeFMDconfig") ||
        (ConfigFn == "DevonSellkeFMDconfig_alternate_transmiss_params")
    InputFilesPrefix = "../../../results/GB_county_model_simn_outputs/Devon_EpiOutputs_Aggregated/"
else
    error("Invalid ConfigFn input.")
end

if (ConfigFn == "CumbriaSellkeFMDconfig") || (ConfigFn == "DevonSellkeFMDconfig")
    BatchID_offset = 15000
elseif (ConfigFn == "CumbriaSellkeFMDconfig_alternate_transmiss_params") ||
        (ConfigFn == "DevonSellkeFMDconfig_alternate_transmiss_params")
    BatchID_offset = 20000
end


#-------------------------------------------------------------------------------
### ITERATE OVER EACH THRESHOLD VALUE TESTED. COMPUTE COST SUMMARY STATISTICS UNDER EACH COST RATIO SCENARIO
#-------------------------------------------------------------------------------

# Specify total number of threshold values to be tested
RiskThresholdVals = collect(0:1:10)
RiskThresholdValsTested = length(RiskThresholdVals)

# Stratification of group proportions
permitted_grp_propn_vals = collect(0:0.05:1)
n_permitted_grp_propn_vals = length(permitted_grp_propn_vals)

# Get total number of scenarios
n_scens = 231

# Declare batch file values to be accessed
# and set an identifier to add to all output file names
for scen_itr = 1:n_scens

    println("scen_itr: $scen_itr")

    # Set range of BatchIDs to be processed
    BatchID_start_idx = BatchID_offset + (scen_itr-1)*RiskThresholdValsTested + 1
    BatchID_end_idx = BatchID_offset + (scen_itr)*RiskThresholdValsTested

    # Set scenario ID
    scen_ID = scen_offset + scen_itr

    # Specify BatchID to be checked
    BatchIDVec = collect(BatchID_start_idx:1:BatchID_end_idx)
    if (ConfigFn == "CumbriaSellkeFMDconfig")
        FileID = "Cumbria_vacc_distance_risk_measure_scenID$scen_ID"
    elseif (ConfigFn == "CumbriaSellkeFMDconfig_alternate_transmiss_params")
        FileID = "Cumbria_alt_pathogen_vacc_distance_risk_measure_scenID$scen_ID"
    elseif (ConfigFn == "DevonSellkeFMDconfig")
        FileID = "Devon_vacc_distance_risk_measure_scenID$scen_ID"
    elseif (ConfigFn == "DevonSellkeFMDconfig_alternate_transmiss_params")
        FileID = "Devon_alt_pathogen_vacc_distance_risk_measure_scenID$scen_ID"
    else
        error("Invalid ConfigFn input.")
    end

    #Specify quantiles to be output (as probability on the interval [0,1])
    QuantileVals = [0.025,0.5,0.975,0,1,0.25,0.75,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9]
        #95% prediction interval, Min, Max, Bounds for 50% prediction interval.
        #Then remaining deciles
    QuantileValsUsed = length(QuantileVals) #Number of quantile values computed

    SummaryStatVals = 1 + length(QuantileVals) #First entry will be the mean

    #Initialise arrays to store summary statistic cost outputs. Row by cost ratio value, column by risk threshold value
    SummStatTotalCost_PremLevel_PopnPersp = zeros(VaccToInfCostRatio_NumTested,RiskThresholdValsTested,SummaryStatVals)
    SummStatTotalCost_AnimalLevel_PopnPersp = zeros(VaccToInfCostRatio_NumTested,RiskThresholdValsTested,SummaryStatVals)

    SummStatTotalCost_PremLevel_IndivPersp = zeros(VaccToInfCostRatio_NumTested,RiskThresholdValsTested,SummaryStatVals)
    SummStatTotalCost_AnimalLevel_IndivPersp = zeros(VaccToInfCostRatio_NumTested,RiskThresholdValsTested,SummaryStatVals)

    #Initialise tuples to store total cost per replicate, under each cost ratio value
    #Tuple entry per risk threshold value
    TotalCostPerSimnReplicate_PremLevel_PopnPersp = Array{Array{Float64,2},1}(undef,RiskThresholdValsTested)
    TotalCostPerSimnReplicate_AnimalLevel_PopnPersp = Array{Array{Float64,2},1}(undef,RiskThresholdValsTested)

    TotalCostPerSimnReplicate_PremLevel_IndivPersp = Array{Array{Float64,2},1}(undef,RiskThresholdValsTested)
    TotalCostPerSimnReplicate_AnimalLevel_IndivPersp = Array{Array{Float64,2},1}(undef,RiskThresholdValsTested)

    #Perform iteration over each threshold value
    for ThresholdIdx = 1:RiskThresholdValsTested

        #---------------------------------------------------------------------------
        ### IMPORT SIMULATION DATA FOR GIVEN THRESHOLD VALUE
        #---------------------------------------------------------------------------
        PremPerDiseaseStateData_EndOfOutbreak = readdlm(string(InputFilesPrefix,"PremPerDiseaseState_BatchID$(BatchIDVec[ThresholdIdx]).txt"))
            # Column order: [S E I R R2 Culled VaccPrem CumulCasePrem CumulCasePremWithNoControlsApplied CumulCasePremControlsApplied]

        CumulInfAnimalsByType = readdlm(string(InputFilesPrefix,"CumulativeCasesAnimals_BatchID$(BatchIDVec[ThresholdIdx]).txt"))
        CumulVaccAnimalsByType = readdlm(string(InputFilesPrefix,"CumulativeVacc_BatchID$(BatchIDVec[ThresholdIdx]).txt"))

        EstIncorrectVaccAnimalsByType = readdlm(string(InputFilesPrefix,"EstIncorrectVaccAnimals_BatchID$(BatchIDVec[ThresholdIdx]).txt"))

        NumAnimalsBothVaccAndInfecByType = readdlm(string(InputFilesPrefix,"NumAnimalsBothVaccAndInfec_BatchID$(BatchIDVec[ThresholdIdx]).txt"))

        #---------------------------------------------------------------------------
        ### ASSIGN IMPORTED DATA TO BE USED IN CALCULATION TO VARIABLES
        #---------------------------------------------------------------------------

        #In all variables, each row corresponds to an individual simulation replicate
        # PremPerDiseaseStateData_EndOfOutbreak column order:
        # [S E I R R2 Culled VaccPrem CumulCasePrem CumulCasePremWithNoControlsApplied CumulCasePremControlsApplied]

        #Infection data, premises level
        PremCumulInfCount = PremPerDiseaseStateData_EndOfOutbreak[:,8]
        PremCumulInfNoControlAppliedCount = PremPerDiseaseStateData_EndOfOutbreak[:,9]

        #Infection data, animal level
        AnimalsCumulInfCount = sum(CumulInfAnimalsByType[:,1:2],dims=2)
            # For each simn replicate (row), sum across each species type the premises
            # that were infected (columns 1 and 2)
        AnimalsCumulInfNoControlAppliedCount = sum(CumulInfAnimalsByType[:,3:4],dims=2)
            # For each simn replicate (row), sum across each species type the premises
            # that were infected and did not apply controls at any stage (columns 3 and 4)

        #Vaccination data
        NumPremVacc = PremPerDiseaseStateData_EndOfOutbreak[:,7]
        NumAnimalsVacc = sum(CumulVaccAnimalsByType,dims=2) #For each simn replicate (row), sum across each species type (so across all columns)

        #Data on premises vaccinated & probability of not being infected over entire outbreak
        NumPremVaccButNoInfec = readdlm(string(InputFilesPrefix,"EstIncorrectVaccPrem_BatchID$(BatchIDVec[ThresholdIdx]).txt"))
        NumAnimalsVaccButNoInfec = sum(EstIncorrectVaccAnimalsByType,dims=2) #For each simn replicate (row), sum across each species type (so across all columns)

        #Data on premises that were vaccinated AND infected
        NumPremBothVaccAndInfec = readdlm(string(InputFilesPrefix,"NumPremBothVaccAndInfec_BatchID$(BatchIDVec[ThresholdIdx]).txt"))
        NumAnimalsBothVaccAndInfec = sum(NumAnimalsBothVaccAndInfecByType,dims=2) #For each simn replicate (row), sum across each species type (so across all columns)

        #---------------------------------------------------------------------------
        ### FOR EACH SIMN REPLICATE, COMPUTE COSTS FOR EACH COST RATIO BY ITERATING OVER EACH COST RATIO
        ### (WITH DATA CORRSPONDING TO A SPECIFIED THRESHOLD VALUE)
        #---------------------------------------------------------------------------

        #Get number of replicates (will match number of rows of infection case count)
        ReplicateNum = length(PremCumulInfCount)

        #Initialise arrays to store total cost per replicate, under each cost ratio value
        #Index of storage tuple given by ThresholdIdx. Tuple entry per risk threshold value tested
        TotalCostPerSimnReplicate_PremLevel_PopnPersp[ThresholdIdx] = zeros(ReplicateNum,VaccToInfCostRatio_NumTested)
        TotalCostPerSimnReplicate_PremLevel_IndivPersp[ThresholdIdx] = zeros(ReplicateNum,VaccToInfCostRatio_NumTested)

        TotalCostPerSimnReplicate_AnimalLevel_PopnPersp[ThresholdIdx] = zeros(ReplicateNum,VaccToInfCostRatio_NumTested)
        TotalCostPerSimnReplicate_AnimalLevel_IndivPersp[ThresholdIdx] = zeros(ReplicateNum,VaccToInfCostRatio_NumTested)

        #Iterate over each cost ratio.
        for CostRatioIdx = 1:VaccToInfCostRatio_NumTested

            #Get the vaccination to infection cost ratio for this iteration
            SelectedVaccToInfCostRatio = VaccToInfCostRatio[CostRatioIdx]

            #National perspective
            TotalCostPerSimnReplicate_PremLevel_PopnPersp[ThresholdIdx][:,CostRatioIdx] = (SelectedVaccToInfCostRatio.*NumPremVacc) .+ PremCumulInfCount # (C_VxNumPremVacc) + (C_IxNumPremInf)
            TotalCostPerSimnReplicate_AnimalLevel_PopnPersp[ThresholdIdx][:,CostRatioIdx] = (SelectedVaccToInfCostRatio.*NumAnimalsVacc) .+ AnimalsCumulInfCount  #(C_VxNumAnimalsVacc) + (C_IxNumAnimalsInf)

            #Farmer perspective
            TotalCostPerSimnReplicate_PremLevel_IndivPersp[ThresholdIdx][:,CostRatioIdx] = (SelectedVaccToInfCostRatio.*NumPremVaccButNoInfec) .+
                                                                                                (1-SelectedVaccToInfCostRatio).*PremCumulInfNoControlAppliedCount .+
                                                                                                (1. *NumPremBothVaccAndInfec) # (C_VxNumPremVaccButNoInfec) + ((C_I-C_V)xNumPremInfAndNoControlApplied) + (C_IxNumPremBothVaccAndInf)
            TotalCostPerSimnReplicate_AnimalLevel_IndivPersp[ThresholdIdx][:,CostRatioIdx] = (SelectedVaccToInfCostRatio.*NumAnimalsVaccButNoInfec) .+
                                                                                                (1-SelectedVaccToInfCostRatio).*AnimalsCumulInfNoControlAppliedCount .+
                                                                                                (1. *NumAnimalsBothVaccAndInfec) #(C_VxNumAnimalsVaccButNoInfec) + ((C_I-C_V)xNumAnimalsInfAndNoControlApplied) + (C_IxNumPremBothVaccAndInf)
        end

        #---------------------------------------------------------------------------
        #  FOR EACH RATIO (AND SINGLE THRESHOLD), TAKE DESIRED COST MEASURE ACROSS SIMULATION REPLICATES (E.G. MEAN, MEDIAN, PERCENTILE,...)
        #---------------------------------------------------------------------------

        #TotalCostPerSimnReplicate arrays are row per simn replicate, column for cost ratio
        #Compute summary statistic by each column

        #First compute the average

        #National perspective
        SummStatTotalCost_PremLevel_PopnPersp[:,ThresholdIdx,1] = mean(TotalCostPerSimnReplicate_PremLevel_PopnPersp[ThresholdIdx],dims=1)
        SummStatTotalCost_AnimalLevel_PopnPersp[:,ThresholdIdx,1] = mean(TotalCostPerSimnReplicate_AnimalLevel_PopnPersp[ThresholdIdx],dims=1)

        #Farmer perspective
        SummStatTotalCost_PremLevel_IndivPersp[:,ThresholdIdx,1] = mean(TotalCostPerSimnReplicate_PremLevel_IndivPersp[ThresholdIdx],dims=1)
        SummStatTotalCost_AnimalLevel_IndivPersp[:,ThresholdIdx,1] = mean(TotalCostPerSimnReplicate_AnimalLevel_IndivPersp[ThresholdIdx],dims=1)

        #Then compute at user defined quantiles
        for QuantileUsedIdx = 1:QuantileValsUsed

            #Quantile to be computed in current loop
            SelectedQuantileVal = QuantileVals[QuantileUsedIdx]

            #Index for the slice the quantile outputs should be assigned to in SummStatTotalCost arrays
            SummStatTotalCostArray_SliceIdx = QuantileUsedIdx + 1

            for CostRatioIdx = 1:VaccToInfCostRatio_NumTested

                #National perspective
                SummStatTotalCost_PremLevel_PopnPersp[CostRatioIdx,ThresholdIdx,SummStatTotalCostArray_SliceIdx] = quantile(TotalCostPerSimnReplicate_PremLevel_PopnPersp[ThresholdIdx][:,CostRatioIdx],SelectedQuantileVal)
                SummStatTotalCost_AnimalLevel_PopnPersp[CostRatioIdx,ThresholdIdx,SummStatTotalCostArray_SliceIdx] = quantile(TotalCostPerSimnReplicate_AnimalLevel_PopnPersp[ThresholdIdx][:,CostRatioIdx],SelectedQuantileVal)

                #Farmer perspective
                SummStatTotalCost_PremLevel_IndivPersp[CostRatioIdx,ThresholdIdx,SummStatTotalCostArray_SliceIdx] = quantile(TotalCostPerSimnReplicate_PremLevel_IndivPersp[ThresholdIdx][:,CostRatioIdx],SelectedQuantileVal)
                SummStatTotalCost_AnimalLevel_IndivPersp[CostRatioIdx,ThresholdIdx,SummStatTotalCostArray_SliceIdx] = quantile(TotalCostPerSimnReplicate_AnimalLevel_IndivPersp[ThresholdIdx][:,CostRatioIdx],SelectedQuantileVal)
            end
        end
    end


    #-------------------------------------------------------------------------------
    ### ITERATE OVER EACH SUMMARY STAISTIC.
    ### FOR EACH RATIO, ACROSS ALL THRESHOLDS CONSIDERED, FIND OPTIMAL THRESHOLD TO MINIMISE COST
    #-------------------------------------------------------------------------------

    #Extract the optimal threshold index (from the cartesian index) and use to get associated risk threshold value
    OptimThresholdVals_PremLevel_PopnPersp = zeros(VaccToInfCostRatio_NumTested,SummaryStatVals)
    OptimThresholdVals_PremLevel_IndivPersp = zeros(VaccToInfCostRatio_NumTested,SummaryStatVals)
    OptimThresholdVals_AnimalLevel_PopnPersp = zeros(VaccToInfCostRatio_NumTested,SummaryStatVals)
    OptimThresholdVals_AnimalLevel_IndivPersp = zeros(VaccToInfCostRatio_NumTested,SummaryStatVals)

    #Iterate over each summary statistic
    for OptimSummStatIdx = 1:SummaryStatVals

        #Get optimal risk threshold values and cartesian indexes (per cost ratio)
        OptimThresholdData_PremLevel_PopnPersp = findmin(SummStatTotalCost_PremLevel_PopnPersp[:,:,OptimSummStatIdx],dims=2)
        OptimThresholdData_PremLevel_IndivPersp = findmin(SummStatTotalCost_PremLevel_IndivPersp[:,:,OptimSummStatIdx],dims=2)
        OptimThresholdData_AnimalLevel_PopnPersp = findmin(SummStatTotalCost_AnimalLevel_PopnPersp[:,:,OptimSummStatIdx],dims=2)
        OptimThresholdData_AnimalLevel_IndivPersp = findmin(SummStatTotalCost_AnimalLevel_IndivPersp[:,:,OptimSummStatIdx],dims=2)

            # Example usage of findmin
            # A = [1.0 2; 3 4]
            #   2×2 Matrix{Float64}:
            #    1.0  2.0
            #    3.0  4.0
            #
            #
            #   julia> findmin(A, dims=2) # Find minimum in each row
            #   ([1.0; 3.0], CartesianIndex{2}[CartesianIndex(1, 1); CartesianIndex(2, 1)])

            # Two tuples. Second tuple is cartesian index

        #Iterate over each cost ratio.
        #Extract the optimal threshold index (from the cartesian index) and use to get assoicated risk threshold value
        for CostRatioIdx = 1:VaccToInfCostRatio_NumTested

            #Each TotalCost array is 2D. Row per cost ratio value, column per risk threshold value.
            #Want the risk threshold value that minimises the total cost at each cost ratio.
            #Therefore, from catesian index tuple, accessed by OptimThresholdData_XXX_XXX[2],
            #find the column index per row (CostRatioIdx) returning the minimum value (Cartesian Index [2])
            # i.e. from cartesian index tuple, access entry [CostRatioIdx][2],
            # giving overall exprssion of OptimThresholdData_XXX_XXX[2][CostRatioIdx][2]
            OptimThresholdIdx_PremLevel_PopnPersp = OptimThresholdData_PremLevel_PopnPersp[2][CostRatioIdx][2]
            OptimThresholdVals_PremLevel_PopnPersp[CostRatioIdx,OptimSummStatIdx] = RiskThresholdVals[OptimThresholdIdx_PremLevel_PopnPersp]

            OptimThresholdIdx_PremLevel_IndivPersp = OptimThresholdData_PremLevel_IndivPersp[2][CostRatioIdx][2]
            OptimThresholdVals_PremLevel_IndivPersp[CostRatioIdx,OptimSummStatIdx] = RiskThresholdVals[OptimThresholdIdx_PremLevel_IndivPersp]

            OptimThresholdIdx_AnimalLevel_PopnPersp = OptimThresholdData_AnimalLevel_PopnPersp[2][CostRatioIdx][2]
            OptimThresholdVals_AnimalLevel_PopnPersp[CostRatioIdx,OptimSummStatIdx] = RiskThresholdVals[OptimThresholdIdx_AnimalLevel_PopnPersp]

            OptimThresholdIdx_AnimalLevel_IndivPersp = OptimThresholdData_AnimalLevel_IndivPersp[2][CostRatioIdx][2]
            OptimThresholdVals_AnimalLevel_IndivPersp[CostRatioIdx,OptimSummStatIdx] = RiskThresholdVals[OptimThresholdIdx_AnimalLevel_IndivPersp]
        end
    end


    #-------------------------------------------------------------------------------
    ### OUTPUT ALL REPLICATE COSTS DATA TO MAT FILE
    #-------------------------------------------------------------------------------

    #Note, in output array:
    # - Row by cost ratio value
    # - Column by risk threshold value

    #Output arrays to file
    AllRepsTotalCostFile = string("../../../results/generate_plot_scripts/OptimBehaviour_AllRepsTotalCostsOutputs/OptimBehav_AllRepsTotalCostData_",FileID,".mat")
    matwrite(AllRepsTotalCostFile, Dict(
            "TotalCostPerSimnReplicate_PremLevel_PopnPersp" => TotalCostPerSimnReplicate_PremLevel_PopnPersp,
            "TotalCostPerSimnReplicate_AnimalLevel_PopnPersp" => TotalCostPerSimnReplicate_AnimalLevel_PopnPersp,
            "TotalCostPerSimnReplicate_PremLevel_IndivPersp" => TotalCostPerSimnReplicate_PremLevel_IndivPersp,
            "TotalCostPerSimnReplicate_AnimalLevel_IndivPersp" => TotalCostPerSimnReplicate_AnimalLevel_IndivPersp,
            ); compress = true)

    #-------------------------------------------------------------------------------
    ### OUTPUT SUMMARY STATISTIC COSTS DATA TO MAT FILE
    #-------------------------------------------------------------------------------

    #Note, in output array:
    # - Row by cost ratio value
    # - Column by risk threshold value

    #Output arrays to file
    SummStatTotalCostFile = string("../../../results/generate_plot_scripts/OptimBehaviour_SummStatsTotalCostsOutputs/OptimBehav_SummStatTotalCostData_",FileID,".mat")
    matwrite(SummStatTotalCostFile, Dict(
            "SummStatTotalCost_PremLevel_PopnPersp" => SummStatTotalCost_PremLevel_PopnPersp,
            "SummStatTotalCost_AnimalLevel_PopnPersp" => SummStatTotalCost_AnimalLevel_PopnPersp,
            "SummStatTotalCost_PremLevel_IndivPersp" => SummStatTotalCost_PremLevel_IndivPersp,
            "SummStatTotalCost_AnimalLevel_IndivPersp" => SummStatTotalCost_AnimalLevel_IndivPersp,
            ); compress = true)

    #-------------------------------------------------------------------------------
    ### OUTPUT THRESHOLD VALUES TO TEXT FILE
    #-------------------------------------------------------------------------------

    # Note, in output array:
    # - Row by cost ratio value
    # - Column by summary statistic

    # Population perspective
    writedlm(string("../../../results/generate_plot_scripts/OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_PremLevel_PopnPersp_",FileID,".txt"),OptimThresholdVals_PremLevel_PopnPersp)
    writedlm(string("../../../results/generate_plot_scripts/OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_PopnPersp_",FileID,".txt"),OptimThresholdVals_AnimalLevel_PopnPersp)

    # Individual perspective
    writedlm(string("../../../results/generate_plot_scripts/OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_PremLevel_IndivPersp_",FileID,".txt"),OptimThresholdVals_PremLevel_IndivPersp)
    writedlm(string("../../../results/generate_plot_scripts/OptimBehaviour_RiskThresholdOutputs/OptimThresholdVals_AnimalLevel_IndivPersp_",FileID,".txt"),OptimThresholdVals_AnimalLevel_IndivPersp)
end
