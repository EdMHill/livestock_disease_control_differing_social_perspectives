#=
Purpose:
File to house functions to carry out local spread algorithms
Used in spatial disease outbreak simulations using the Sellke construction

Algorithm list:
- Pairwise algorithm (single susceptible unit targeted by a single infectious unit)

Date: 3rd November 2021
=#

"""
    ComputeLocalPairwiseProb(SusNode_xLoc::Float64,
                            SusNode_yLoc::Float64,
                            InfNode_xLoc::Float64,
                            InfNode_yLoc::Float64,
                            PremTransmiss::Array{Float64},
                            PremSuscept::Array{Float64},
                            CoordType::Int64,
                            KernelLookUpVec::Array{Float64,1},
                            delta_t::Float64,
                            SelectedInfNodeID::Int64,
                            SelectedSusNodeID::Int64)

Calculates directly for each infectious-susceptible pair in the population the probability of infection occurring, which is evaluated as a Benoulli process.

Inputs:
- `SusNode_xLoc::Array{Float64,1}`: East-west plane co-ordinate of susceptible premises.
- `SusNode_yLoc::Array{Float64,1}`: North-south plane co-ordinate of susceptible premises.
- `InfNode_xLoc::Array{Float64,1}`: East-west plane co-ordinate of infectious premises.
- `InfNode_yLoc::Array{Float64,1}`: North-south plane co-ordinate of infectious premises.
- `PremTransmiss::Array{Float64}`: Premises-level transmissibility.
- `PremSuscept::Array{Float64}`: Premises-level susceptibility.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `KernelLookUpVec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `delta_t::Float64`: Timestep increment.
- `SelectedInfNodeID::Int64`: ID of the infectious node.
- `SelectedSusNodeID::Int64`: ID of the susceptible node.

Outputs:
- `PairwiseInfProb::Float64`: Probability of infection.

Location: SellkeLocalSpreadFns.jl
"""
function ComputeLocalPairwiseProb(SusNode_xLoc::Float64,
                                    SusNode_yLoc::Float64,
                                    InfNode_xLoc::Float64,
                                    InfNode_yLoc::Float64,
                                    PremTransmiss::Array{Float64},
                                    PremSuscept::Array{Float64},
                                    CoordType::Int64,
                                    KernelLookUpVec::Array{Float64,1},
                                    delta_t::Float64,
                                    SelectedInfNodeID::Int64,
                                    SelectedSusNodeID::Int64)

    #Calculate distance between the two points
    if CoordType == 1 #Cartesian co-ords (metres)
        d =  eucl_distance(InfNode_xLoc,
                            InfNode_yLoc,
                            SusNode_xLoc,
                            SusNode_yLoc)
    elseif CoordType == 2 #Cartesian co-ords (metres)
        d = eucl_distance_ConvertToMetres(InfNode_xLoc,
                                            InfNode_yLoc,
                                            SusNode_xLoc,
                                            SusNode_yLoc)
    elseif CoordType == 3 #Lat/Long co-ords
        println("In LatLong loop!")
        d = GreatCircleDistance(InfNode_yLoc, InfNode_xLoc,  #lat1, lon1
                                            SusNode_yLoc, SusNode_xLoc) #lat2, lon2
    end

    #Calculate rate of infection
    distIdx = ReturnDistIdxForKernel(d)
    FOIrate = PremTransmiss[SelectedInfNodeID]*PremSuscept[SelectedSusNodeID]*KernelLookUpVec[distIdx]

    #Compute probability of infection (Use 1-exp function!)
    PairwiseInfProb = oneMinusExp(-FOIrate*delta_t)

    return PairwiseInfProb::Float64
end
