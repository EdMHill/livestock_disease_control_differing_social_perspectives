#=
Purpose:
File to house functions to compute distances between point locations
Used in spatial disease outbreak simulations

Formats:
 - Euclidean distance squared
 - Euclidean distance
 - Great circle distance

Date: 3rd November 2021
=#

#-------------------------------------------------------------------------------
### EUCLIDEAN DISTANCE SQUARED (BETWEEN NODES)
#-------------------------------------------------------------------------------
"""
    sq_distance(Loc_A_xVal::Float64,
                Loc_A_yVal::Float64,
                Loc_B_xVal::Float64,
                Loc_B_yVal::Float64)

Return squared euclidean distance between two nodes, in metres.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: DistanceFns.jl
"""
function sq_distance(Loc_A_xVal::Float64,
                            Loc_A_yVal::Float64,
                            Loc_B_xVal::Float64,
                            Loc_B_yVal::Float64)

    delta_x::Float64 = Loc_A_xVal - Loc_B_xVal #Find difference in east-west plane
    delta_y::Float64 = Loc_A_yVal - Loc_B_yVal #Find difference in south-north plane

    return (delta_x*delta_x + delta_y*delta_y)::Float64
end

"""
    sq_distance_ConvertToKM(Loc_A_xVal::Float64,
                            Loc_A_yVal::Float64,
                            Loc_B_xVal::Float64,
                            Loc_B_yVal::Float64)

Return squared euclidean distance between two nodes, in km.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: DistanceFns.jl
"""
function sq_distance_ConvertToKM(Loc_A_xVal::Float64,
                            Loc_A_yVal::Float64,
                            Loc_B_xVal::Float64,
                            Loc_B_yVal::Float64)

    return ((sq_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal))*0.001)::Float64
end

"""
    sq_distance_ConvertToMetres(Loc_A_xVal::Float64,
                                Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64,
                                Loc_B_yVal::Float64)

Return squared euclidean distance between two nodes, in metres.\n
Used when inputs using a co-ordinate system with distances in units of km.

Location: DistanceFns.jl
"""
function sq_distance_ConvertToMetres(Loc_A_xVal::Float64,
                            Loc_A_yVal::Float64,
                            Loc_B_xVal::Float64,
                            Loc_B_yVal::Float64)

    return ((sq_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal))*1000.)::Float64
end

#-------------------------------------------------------------------------------
### EUCLIDEAN DISTANCE (BETWEEN NODES)
#-------------------------------------------------------------------------------
"""
    eucl_distance(Loc_A_xVal::Float64,
                Loc_A_yVal::Float64,
                Loc_B_xVal::Float64,
                Loc_B_yVal::Float64)

Return euclidean distance between two nodes, in metres.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: DistanceFns.jl
"""
function eucl_distance(Loc_A_xVal::Float64, Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64, Loc_B_yVal::Float64)


    # Take square root of squared distance between nodes A & B
    return (sqrt(sq_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal)))::Float64

end

"""
    eucl_distance_ConvertToKM(Loc_A_xVal::Float64,
                                Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64,
                                Loc_B_yVal::Float64)

Return euclidean distance between two nodes, in km.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: DistanceFns.jl
"""
function eucl_distance_ConvertToKM(Loc_A_xVal::Float64, Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64, Loc_B_yVal::Float64)


    return ((eucl_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal))*0.001)::Float64

end

"""
    eucl_distance_ConvertToMetres(Loc_A_xVal::Float64,
                                    Loc_A_yVal::Float64,
                                    Loc_B_xVal::Float64,
                                    Loc_B_yVal::Float64)

Return euclidean distance between two nodes, in metres.\n
Used when inputs using a co-ordinate system with distances in units of km.

Location: DistanceFns.jl
"""
function eucl_distance_ConvertToMetres(Loc_A_xVal::Float64, Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64, Loc_B_yVal::Float64)


    return ((eucl_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal))*1000.)::Float64

end

#-------------------------------------------------------------------------------
### GREAT CIRCLE DISTANCE (BY DEFAULT IN KM, SCALE UP TO METRES)
#-------------------------------------------------------------------------------

"""
    GreatCircleDistance(lat1::Float64, lon1::Float64, lat2::Float64, lon2::Float64)

Use the haversine formula to give the great-circle distances between two points on a sphere from their longitudes and latitudes.

Note, unit of output is metres.

Location: DistanceFns.jl
"""
function GreatCircleDistance(lat1::Float64, lon1::Float64, lat2::Float64, lon2::Float64)

    #Compute haversine formula
    LatLongDist = 2 * 6371 * asin(sqrt(sind((lat2 - lat1)*0.5)*sind((lat2 - lat1)*0.5) +
                                            cosd(lat1) * cosd(lat2) * sind((lon2 - lon1)*0.5)* sind((lon2 - lon1)*0.5)))
    #LatLongDist = 2 * 6371 * asin(sqrt(sind((lat2 - lat1)*0.5)^2 +
    #                                cosd(lat1) * cosd(lat2) * sind((lon2 - lon1)*0.5)^2))

    #Return distance
    return LatLongDist*1000.::Float64

end

"""
    GreatCircleDistanceBetweenCells(Cell_A_xMin::Float64, Cell_A_xMax::Float64,
                                    Cell_A_yMin::Float64, Cell_A_yMax::Float64,
                                    Cell_B_xMin::Float64, Cell_B_xMax::Float64,
                                    Cell_B_yMin::Float64, Cell_B_yMax::Float64)

Use the haversine formula to give the shortest great-circle distances between two cells on a sphere from their longitudes and latitudes.

Note, unit of output is metres.

Location: DistanceFns.jl
"""
function GreatCircleDistanceBetweenCells(Cell_A_xMin::Float64, Cell_A_xMax::Float64,
                                            Cell_A_yMin::Float64, Cell_A_yMax::Float64,
                                            Cell_B_xMin::Float64, Cell_B_xMax::Float64,
                                            Cell_B_yMin::Float64, Cell_B_yMax::Float64)

    #Get shortest horizontal distance between cells A & B
    #Ensure all distances are positive, use absolute value fn
    CellDistances_x = abs.([Cell_A_xMin - Cell_B_xMin,  #Cell_A.xVal_min - Cell_B.xVal_min
                            Cell_A_xMin - Cell_B_xMax,  #Cell_A.xVal_min - Cell_B.xVal_max
                            Cell_A_xMax - Cell_B_xMin,  #Cell_A.xVal_max - Cell_B.xVal_min
                            Cell_A_xMax - Cell_B_xMax])  #Cell_A.xVal_max - Cell_B.xVal_max

    MinGridDist_xIdx::Int64 = findmin(CellDistances_x)[2] #Second entry, index of the minimum over the given dimensions.


    #Will take value 1, 2, 3 or 4. Throw error if another value
    #Lookup table based on index value obtained
    #Assign longitude coords
    if MinGridDist_xIdx == 1 #Cell_A min, Cell_B min
        lon1 = Cell_A_xMin
        lon2 = Cell_B_xMin
    elseif MinGridDist_xIdx == 2 #Cell_A min, Cell_B max
        lon1 = Cell_A_xMin
        lon2 = Cell_B_xMax
    elseif MinGridDist_xIdx == 3 #Cell_A max, Cell_B min
        lon1 = Cell_A_xMax
        lon2 = Cell_B_xMin
    elseif MinGridDist_xIdx == 4 #Cell_A max, Cell_B max
        lon1 = Cell_A_xMax
        lon2 = Cell_B_xMax
    else #Unexpected value for MinGridDist_xIdx, throw error
        error("MinGridDist_xIdx has value $MinGridDist_xIdx. MinGridDist_xIdx must have value 1, 2, 3 or 4.")
    end

    #Get shortest vertical distance between cells A & B
    #Ensure all distances are positive, use absolute value fn
    CellDistances_y = abs.([Cell_A_yMin - Cell_B_yMin,  #Cell_A.yVal_min - Cell_B.yVal_min
                            Cell_A_yMin - Cell_B_yMax,  #Cell_A.yVal_min - Cell_B.yVal_max
                            Cell_A_yMax - Cell_B_yMin,  #Cell_A.yVal_max - Cell_B.yVal_min
                            Cell_A_yMax - Cell_B_yMax])  #Cell_A.yVal_max - Cell_B.yVal_max

    MinGridDist_yIdx::Int64 = findmin(CellDistances_y)[2] #Second entry, index of the minimum over the given dimensions.


    #Will take value 1, 2, 3 or 4. Throw error if another value
    #Lookup table based on index value obtained
    #Assign latitude coords
    if MinGridDist_yIdx == 1 #Cell_A min, Cell_B min
        lat1 = Cell_A_yMin
        lat2 = Cell_B_yMin
    elseif MinGridDist_yIdx == 2 #Cell_A min, Cell_B max
        lat1 = Cell_A_yMin
        lat2 = Cell_B_yMax
    elseif MinGridDist_yIdx == 3 #Cell_A max, Cell_B min
        lat1 = Cell_A_yMax
        lat2 = Cell_B_yMin
    elseif MinGridDist_yIdx == 4 #Cell_A max, Cell_B max
        lat1 = Cell_A_yMax
        lat2 = Cell_B_yMax
    else #Unexpected value for MinGridDist_xIdx, throw error
        error("MinGridDist_yIdx has value $MinGridDist_yIdx. MinGridDist_yIdx must have value 1, 2, 3 or 4.")
    end

    #Return great circle distance
    return (GreatCircleDistance(lat1, lon1, lat2, lon2))::Float64

end
