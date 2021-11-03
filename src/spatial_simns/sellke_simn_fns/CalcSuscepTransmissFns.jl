#=
Purpose:
File to house functions used for calculating premises-level susceptibility and transmissibility

Date: 3rd November 2021
=#

#-------------------------------------------------------------------------------
### SINGLE ANIMAL TYPE
#-------------------------------------------------------------------------------

"""
    CalcPremSuscepTransmiss_FMDlike(PremLivestockData::Array{Int64},
                                        PerAnimalSuscep::Array{Float64,1},
                                        SuscepExponent::Array{Float64,1},
                                        PerAnimalTransmiss::Array{Float64,1},
                                        TransmissExponent::Array{Float64,1})

Calculate premises-level and livestock type level susceptibility and transmissibility for pathogen with characteristics akin to FMD.

Inputs:
- `PremLivestockData::Array{Int64}`: Number of each livestock type per premises. Row per premises, column per animal.
- `PerAnimalSuscep::Array{Float64,1}`: Susceptibility scale parameters for each species.
- `SuscepExponent::Array{Float64,1}`: Transmissibility scale parameters for each species.
- `PerAnimalTransmiss::Array{Float64,1}`: Susceptibility exponents for each species.
- `TransmissExponent::Array{Float64,1}`: Transmissibility exponents for each species.

Outputs:
- `PremSuscept::Array{Float64,1}`, `PremTransmiss::Array{Float64,1}`: Premises-level susceptibility and transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`, `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, livestock type level susceptibility and transmissibility

Location: CalcSuscepTransmissFns.jl
"""
function CalcPremSuscepTransmiss_FMDlike(PremLivestockData::Array{Int64},
                                            PerAnimalSuscep::Array{Float64,1}, SuscepExponent::Array{Float64,1},
                                            PerAnimalTransmiss::Array{Float64,1}, TransmissExponent::Array{Float64,1})

    #Use PremLivestockData input data,(Number of each livestock type per premises. Row per premises, column per animal)
    #with PerAnimalSuscep, PerAnimalSuscep, SuscepExponent, TransmissExponent
    #For each premises, get contribution to susceptibility and transmissibility from each species
    PremSuscept_SpeciesBreakdown = PerAnimalSuscep.*(PremLivestockData.^SuscepExponent)
    PremTransmiss_SpeciesBreakdown = PerAnimalTransmiss.*(PremLivestockData.^TransmissExponent)

    #Sum across columns to get overall premises-level value
    PremSuscept = sum(PremSuscept_SpeciesBreakdown,dims=2)
    PremTransmiss = sum(PremTransmiss_SpeciesBreakdown,dims=2)

    #Calculation outline:
    #Exponent i applied to col_i of PremLivestockData
    #Multiply col_i of PerAnimalSuscep/PerAnimalTransmiss by col_i of PremLivestockData
    #Sum across columns to get overall premises-level value

    return PremSuscept::Array{Float64,1}, PremTransmiss::Array{Float64,1},
            PremSuscept_SpeciesBreakdown::Array{Float64,2}, PremTransmiss_SpeciesBreakdown::Array{Float64,2}
end

"""
    CalcPremSuscepTransmiss_FMDlike(PremLivestockData::Array{Int64},
                                        PerAnimalSuscep::Array{Float64,2},
                                        SuscepExponent::Array{Float64,2},
                                        PerAnimalTransmiss::Array{Float64,2},
                                        TransmissExponent::Array{Float64,2})

Calculate premises-level and livestock type level susceptibility and transmissibility for pathogen with characteristics akin to FMD.

Inputs:
- `PremLivestockData::Array{Int64}`: Number of each livestock type per premises. Row per premises, column per animal.
- `PerAnimalSuscep::Array{Float64,1}`: Susceptibility scale parameters for each species.
- `SuscepExponent::Array{Float64,1}`: Transmissibility scale parameters for each species.
- `PerAnimalTransmiss::Array{Float64,1}`: Susceptibility exponents for each species.
- `TransmissExponent::Array{Float64,1}`: Transmissibility exponents for each species.

Outputs:
- `PremSuscept::Array{Float64,2}`, `PremTransmiss::Array{Float64,2}`: Premises-level susceptibility and transmissibility.
- `PremSuscept_SpeciesBreakdown::Array{Float64,2}`, `PremTransmiss_SpeciesBreakdown::Array{Float64,2}`: Per premises, livestock type level susceptibility and transmissibility

Location: CalcSuscepTransmissFns.jl
"""
function CalcPremSuscepTransmiss_FMDlike(PremLivestockData::Array{Int64},
                                            PerAnimalSuscep::Array{Float64,2}, SuscepExponent::Array{Float64,2},
                                            PerAnimalTransmiss::Array{Float64,2}, TransmissExponent::Array{Float64,2})

    #Use PremLivestockData input data,(Number of each livestock type per premises. Row per premises, column per animal)
    #with PerAnimalSuscep, PerAnimalSuscep, SuscepExponent, TransmissExponent
    #For each premises, get contribution to susceptibility and transmissibility from each species
    PremSuscept_SpeciesBreakdown = PerAnimalSuscep.*(PremLivestockData.^SuscepExponent)
    PremTransmiss_SpeciesBreakdown = PerAnimalTransmiss.*(PremLivestockData.^TransmissExponent)

    #Sum across columns to get overall premises-level value
    PremSuscept = sum(PremSuscept_SpeciesBreakdown,dims=2)
    PremTransmiss = sum(PremTransmiss_SpeciesBreakdown,dims=2)

    #Calculation outline:
    #Exponent i applied to col_i of PremLivestockData
    #Multiply col_i of PerAnimalSuscep/PerAnimalTransmiss by col_i of PremLivestockData
    #Sum across columns to get overall premises-level value

    return PremSuscept::Array{Float64,2}, PremTransmiss::Array{Float64,2},
            PremSuscept_SpeciesBreakdown::Array{Float64,2}, PremTransmiss_SpeciesBreakdown::Array{Float64,2}
end
