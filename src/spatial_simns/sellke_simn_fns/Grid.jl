#=
Purpose:
File to house functions used to seed infection and check landscape parameters are valid

Alogorithm list:
- Seed infection
- Check landscape & premises data are compatible
- Return landscaspe size attributes (& check values are valid)
- Return distance index, used to look up Kernel value

Date: 3rd November 2021
=#

#-------------------------------------------------------------------------------
### METHODS OF SEEDING INFECTION
#-------------------------------------------------------------------------------
"""
    SeedInfection(InitialInfInfo,
                    PremNum::Int64,
                    PremLoc_AllVals::Array{Float64,2},
                    CoordType::Int64,
                    rng::MersenneTwister)

Seed initial infection cases according to a specified method.

Inputs:
- `InitialInfInfo`: (tuple) [SeedMethod,NumOfNodes/NodeIDs]
                             Seed method.
                              - 1 = random,
                              - 2/3 = single/group specific node id(s),
                              - 4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
                              - 5 = seed a random site and it's N nearest neighbours.
                             Number of nodes to seed each replicate if SeedMethod = 1 or 5; or if SeedMethod = 2/3, seed this specific node every replicate.
                                SeedMethod = 4 from file, node ids to seed given by one id/row, number of lines must be == number of replicates.
- `PremNum::Int64`: Number of premises in landscape.
- `PremLoc_AllVals::Array{Float64,2}`: Coordinates for each premises.
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `rng::AbstractRNG`: The random number generator.

Outputs:
- `SeedPremIDs::Array{Int64,1}`: Vector of 1's (infected) and 0's (not infected). Entry per node that has infection.

Location: Grid.jl
"""
function SeedInfection(InitialInfInfo,
                        PremNum::Int64,
                        PremLoc_AllVals::Array{Float64,2},
                        CoordType::Int64,
                        rng::AbstractRNG)

    #Disaggregate InitialInfInfo
    SeedMethod = InitialInfInfo[1]::Int64


    if SeedMethod == 1 #Random
        #Get number of premises to be seeded with infection each replicate
        NumPremToSeed = InitialInfInfo[2]::Int64

        # Draw IDs based on number of premises in use.
        # Use rng to get consistency in selected premises IDs across control runs
        SeedPremIDs = rand(rng,1:PremNum,NumPremToSeed)
    elseif SeedMethod == 2 #Single specific node ID
        SeedPremIDs = [InitialInfInfo[2]]::Array{Int64,1}
    elseif SeedMethod == 3 #Set of specific node ID
        SeedPremIDs = InitialInfInfo[2]::Array{Int64,1}
    elseif SeedMethod == 4 #From file (filename stored in InitialInfInfo[2])
        #One id per row, will seed all the node ids given in the file each replicate.
        SeedPremIDs = readdlm(InitialInfInfo[2], Int64)
    elseif SeedMethod == 5 #Single random node and it's (NumPremToSeed-1) nearest neighbours

        # Disaggregate location data
        PremLoc_xVals::Array{Float64,1} = PremLoc_AllVals[:,1]
        PremLoc_yVals::Array{Float64,1} = PremLoc_AllVals[:,2]

        #Get number of premises to be seeded with infection each replicate
        NumPremToSeed = InitialInfInfo[2]::Int64

        # Draw the index case from which cluster will be generated
        # Use rng to get consistency in selected premises IDs across control runs
        index_SeedPremID = rand(rng,1:PremNum,1)[1]
        println("index_SeedPremID: $index_SeedPremID")

        # Get distances from index case to all other premises
        # Check distance to premises now reporting infection
        dist_to_index_prem = zeros(Float64,PremNum) # Initialise vector to store each premises to index premises distance
        for Prem_ID = 1:PremNum
            if CoordType == 1 #Cartesian co-ords (metres)
                dist_to_index_prem[Prem_ID] =  eucl_distance(PremLoc_xVals[index_SeedPremID],
                                    PremLoc_yVals[index_SeedPremID],
                                    PremLoc_xVals[Prem_ID],
                                    PremLoc_yVals[Prem_ID])
            elseif CoordType == 2 #Cartesian co-ords (metres)
                dist_to_index_prem[Prem_ID] = eucl_distance_ConvertToMetres(PremLoc_xVals[index_SeedPremID],
                                                    PremLoc_yVals[index_SeedPremID],
                                                    PremLoc_xVals[Prem_ID],
                                                    PremLoc_yVals[Prem_ID])
            elseif CoordType == 3 #Lat/Long co-ords
                println("In LatLong loop!")
                dist_to_index_prem[Prem_ID] = GreatCircleDistance(PremLoc_yVals[index_SeedPremID], PremLoc_xVals[index_SeedPremID],  #lat1, lon1
                                                    PremLoc_yVals[Prem_ID], PremLoc_xVals[Prem_ID]) #lat2, lon2
            end
        end

        # Get indexes of sorted ascending order
        ordered_distance_prem_IDs = sortperm(dist_to_index_prem)

        # Remove the index case ID
        filter!(x->x≠index_SeedPremID,ordered_distance_prem_IDs)

        # Retain IDs of index case & (NumPremToSeed-1) nearest neighbours
        SeedPremIDs = [index_SeedPremID;ordered_distance_prem_IDs[1:(NumPremToSeed-1)]]

        # Error check
        # Should be no duplicates in SeedPremIDs
        if length(unique(SeedPremIDs)) < NumPremToSeed
            error("length(unique(SeedPremIDs)): $(length(unique(SeedPremIDs)))); Less than NumPremToSeed ($NumPremToSeed). Invalid.")
        end
    else #Invalid value, throw error
        error("SeedMethod has value $SeedMethod. SeedMethod should take value 1, 2, 3, 4 or 5.")
    end

    return SeedPremIDs::Array{Int64,1}
end

#-------------------------------------------------------------------------------
### CHECK LANDSCAPE & PREMISES DATA ARE COMPATIBLE
#-------------------------------------------------------------------------------
"""
    CheckLandscapeValid(BoundingBoxVar::Array{Float64,1},
                        PremLoc_xVals::Array{Float64,1},
                        PremLoc_yVals::Array{Float64,1})

Check landscape and premises data are compatible.

Inputs:
- `BoundingBoxVar::Array{Float64,1}`: Vector with limits for landscape bounding box. Entries [Min_x,Max_x,Min_y,Max_y].
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).

Outputs: None \n
Location: Grid.jl
"""
function CheckLandscapeValid(BoundingBoxVar::Array{Float64,1},
                            PremLoc_xVals::Array{Float64,1},
                            PremLoc_yVals::Array{Float64,1})

    #Function designed to throw error and exit programme if premises lie on boundary/outside
    #landscape bounding box

    #CHECK BOUNDING BOX WIDTH AND HEIGHT WILL BE POSITIVE!
    xMin =  BoundingBoxVar[1]; xMax = BoundingBoxVar[2];
    yMin =  BoundingBoxVar[3]; yMax = BoundingBoxVar[4];
    BoundingBoxWidth = xMax - xMin
    BoundingBoxHeight = yMax - yMin

    println("PremLoc_xVals min: $(minimum(PremLoc_xVals))")
    println("PremLoc_xVals max: $(maximum(PremLoc_xVals))")
    println("PremLoc_yVals min: $(minimum(PremLoc_yVals))")
    println("PremLoc_yVals max: $(maximum(PremLoc_yVals))")

    if BoundingBoxWidth <= 0 && BoundingBoxHeight <= 0
        error("Nonpositive BoundingBoxWidth and BoundingBoxHeight found.")
    elseif BoundingBoxWidth <= 0
        error("Nonpositive BoundingBoxWidth found.")
    elseif BoundingBoxHeight <= 0
        error("Nonpositive BoundingBoxHeight found.")
    end

    #CHECK THAT NO PREMISES WILL BE LOCATED ON AN EDGE OR OUTSIDE THE BOUNDING BOX
    #Locations on bottom edge and left edge are okay, such premises will be allocated a grid in the grid configuration.
    #Premises on top edge and right edge would be missed!
    ErrorFlag = 0 #Initialise variable to determine whether error should be thrown or not
    if (sum(PremLoc_xVals .<= xMin) > 0)
        println("Incompatible boundary. Premises located on or outside left edge of landscape box.")
        ErrorFlag = 1
    end

    if (sum(PremLoc_xVals .>= xMax) > 0)
        println("Incompatible boundary. Premises located on or outside right edge of landscape box.")
        ErrorFlag = 1
    end

    if (sum(PremLoc_yVals .<= yMin) > 0)
        println("Incompatible boundary. Premises located on or below bottom edge of landscape box.")
        ErrorFlag = 1
    end

    if (sum(PremLoc_yVals .>= yMax) > 0)
        println("Incompatible boundary. Premises located on or above top edge of landscape box.")
        ErrorFlag = 1
    end

    if ErrorFlag == 1
        error("Incompatible landscape bounding box. Amend boundary values. Exiting programme.")
    end

    return
end
#-------------------------------------------------------------------------------
### GET LANDSCAPE SIZE ATTRIBUTES (& CHECK VALUES ARE VALID)
#-------------------------------------------------------------------------------
"""
    GetLandscapeSize(BoundingBoxVar::Array{Float64,1},
                        PremLoc_xVals::Array{Float64,1},
                        PremLoc_yVals::Array{Float64,1},
                        CoordType::Int64)

Get landscape attributes and check the attributes are valid.

Inputs:
- `BoundingBoxVar::Array{Float64,1}`: Vector with limits for landscape bounding box. Entries [Min_x,Max_x,Min_y,Max_y].
- `PremLoc_xVals::Array{Float64,1}`: East-west plane co-ordinate (per premises).
- `PremLoc_yVals::Array{Float64,1}`: North-south plane co-ordinate (per premises).
- `CoordType::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").

Outputs:
- `BoundingBoxWidth::Float64`: Width of the bounding box for the landscape.
- `BoundingBoxHeight::Float64`: Height of the bounding box for the landscape.
- `LongestLandscapeEdge::Float64`: Longest edge of the bounding box for the landscape.
- `MaxDistWithinLandscape::Int64`: Distance of diagonal across landscape (maximum separation between two points within the landscape).

Location: Grid.jl
"""
function GetLandscapeSize(BoundingBoxVar::Array{Float64,1},
                            PremLoc_xVals::Array{Float64,1},
                            PremLoc_yVals::Array{Float64,1},
                            CoordType::Int64)
    #Calculate length of longest edge
    BoundingBoxWidth = BoundingBoxVar[2] - BoundingBoxVar[1]
    BoundingBoxHeight = BoundingBoxVar[4] - BoundingBoxVar[3]
    LongestLandscapeEdge = max(BoundingBoxWidth,BoundingBoxHeight)
    println("LongestLandscapeEdge is $LongestLandscapeEdge. CHECK THIS IS CORRECT!")

    #Calculate distance (in metres) along diagonal for landscape
    if CoordType == 1 #Cartesian co-ordinate system. Gives distance between nodes in metres.
        MaxDistWithinLandscape_float = sqrt(LongestLandscapeEdge*LongestLandscapeEdge + LongestLandscapeEdge*LongestLandscapeEdge)
    elseif CoordType == 2 #Cartesian co-ordinate system. Gives distance between nodes in km. Convert to metres.
        MaxDistWithinLandscape_km = sqrt(LongestLandscapeEdge*LongestLandscapeEdge + LongestLandscapeEdge*LongestLandscapeEdge)
        MaxDistWithinLandscape_float = MaxDistWithinLandscape_km*1000
    else #LatLong co-ordinate system. Output in km. Convert to metres.
        MaxDistWithinLandscape_float = GreatCircleDistance(BoundingBoxVar[1]+LongestLandscapeEdge, BoundingBoxVar[1],
                                                            BoundingBoxVar[2]+LongestLandscapeEdge, BoundingBoxVar[2])
    end

    #Take ceiling of float distance and convert to integer
    MaxDistWithinLandscape = convert(Int64,ceil(MaxDistWithinLandscape_float))
    println("MaxDistWithinLandscape is $MaxDistWithinLandscape. CHECK THIS IS CORRECT!")

    return BoundingBoxWidth::Float64,
            BoundingBoxHeight::Float64,
            LongestLandscapeEdge::Float64,
            MaxDistWithinLandscape::Int64
end

#-------------------------------------------------------------------------------
### RETURN DISTANCE INDEX, USED TO LOOK UP KERNEL VALUE
#-------------------------------------------------------------------------------
"""
    ReturnDistIdxForKernel(d::Float64)

Return distance index, used to look up kernel value.

Inputs:
- `d::Float64`: Distance between epidemiological units of interest.

Outputs:
- `distIdx::Int64`: Array index to use in kernel lookup array.

Location: Grid.jl
"""
function ReturnDistIdxForKernel(d::Float64)

    if d <= 0.5
        distIdx = 1
    else  #For distances of 0m, use first entry in vector.
        distIdx = round(Int, d)
    end

    return distIdx::Int64
end
