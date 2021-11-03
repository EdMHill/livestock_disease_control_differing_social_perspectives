#=
Purpose:
File to house kernel functions
Used in spatial disease outbreak simulations

Kernel list:
1: Buhnerkempe (from Buhnerkempe et al 2014 & used in Sellman et al 2018)
2: Buhnerkempe half gamma
3: Buhnerkempe double gamma
4: gaussian
5: Brand a3
6: Brand a4
7: Brand a5
8: Hayama (from Hayama et al 2015)
9: USDOSv2 (from Tsao et al 2020)
10: GB FMD kernel (from Keeling & Rohani)

Date: 3rd November 2021
=#

# #Generic format of kernel functions
# #Accept float as input. Return float as output
# function(d::float)
#
#     #List of required parameter values
#     k1 = 0
#     k2 = 0
#     ...
#
#     #Express kernel formulation
#     KernelVal = K(d)
#
#     #Return kernal value for distance d from function
#     return KernelVal
#
# end


#-------------------------------------------------------------------------------
### 1: Buhnerkempe kernel (from Buhnerkempe et al 2014 & used in Sellman et al 2018)
#-------------------------------------------------------------------------------

"""
    Construct_USDOS2_kernel(MaxDist::Int64)

Create lookup vector for Buhnerkempe transmission kernel (from Buhnerkempe et al 2014). Entry for each one metre increment.

Inputs:
- `MaxDist::Int64`: Maximum distance to compute value of transmission kernel for.

Outputs:
- `KernelLookUpVec::Array{Float64,1}`: Kernel term, with entry for each one metre increment.

Location: SpatialKernels.jl
"""
function Construct_USDOS2_kernel(MaxDist::Int64)

    #Initialise lookup array
    KernelLookUpVec = zeros(Float64,MaxDist)

    #Iterate over desired distances. Assign kernel value to array
    for DistItr = 1:MaxDist
        KernelLookUpVec[DistItr] = USDOS2_kernel(convert(Float64,DistItr))
    end

    return KernelLookUpVec::Array{Float64,1}
end

"""
    USDOS2_kernel(d::Float64)

Calculate value of Buhnerkempe transmission kernel (from Buhnerkempe et al 2014).

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function USDOS2_kernel(d::Float64)

    #List of required parameter values
    k1::Float64 = 0.08912676
    k2::Float64 = 1600.0
    k3::Float64 = 4.6

    #Express kernel formulation
    KernelVal::Float64 = k1 / (1 + ((d/k2)^k3))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end

"""
    USDOS2_kernel_square(d::Float64)

Calculate value of Buhnerkempe transmission kernel (from Buhnerkempe et al 2014).\n
Same as USDOS2\\_kernel but accepts distance^2 instead of distance.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function USDOS2_kernel_square(d::Float64)

    # Same as USDOS2_kernel but accepts distance^2 instead of distance.
    # This means that sqrt operation can be avoided when calc. distance,
    # which potentially could be a little bit faster.

    #List of required parameter values
    k1::Float64 = 0.08912676
    k2::Float64 = 1600.0
    k3::Float64 = 4.6

    #Express kernel formulation
    KernelVal::Float64 = k1 / (1 + ((d^(0.5*k3))/(k2^k3)))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end

#-------------------------------------------------------------------------------
### 2: Buhnerkempe half gamma kernels
#-------------------------------------------------------------------------------

"""
    USDOS2_kernel_half_gamma(d::Float64)

Calculate value of Buhnerkempe half gamma transmission kernel.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function USDOS2_kernel_half_gamma(d::Float64)

    #List of required parameter values
    k1::Float64 = 0.01813338
    k2::Float64 = 1600.0
    k3::Float64 = 2.3

    #Express kernel formulation
    KernelVal::Float64 = k1 / (1 + ((d/k2)^k3))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end


"""
    USDOS2_kernel_square_half_gamma(d::Float64)

Calculate value of Buhnerkempe half gamma transmission kernel.\n
Accepts distance^2 instead of distance.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function USDOS2_kernel_square_half_gamma(d::Float64)

    #List of required parameter values
    k1::Float64 = 0.01813338
    k2::Float64 = 1600.0
    k3::Float64 = 2.3

    #Express kernel formulation
    KernelVal::Float64 = k1 / (1 + ((d^(0.5*k3))/(k2^k3)))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end

#-------------------------------------------------------------------------------
### 3: Buhnerkempe double gamma kernels
#-------------------------------------------------------------------------------

"""
    USDOS2_kernel_double_gamma(d::Float64)

Calculate value of Buhnerkempe double gamma transmission kernel.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function USDOS2_kernel_double_gamma(d::Float64)

    #List of required parameter values
    k1::Float64 = 0.11489682
    k2::Float64 = 1600.0
    k3::Float64 = 9.2

    #Express kernel formulation
    KernelVal::Float64 = k1 / (1 + ((d/k2)^k3))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end


"""
    USDOS2_kernel_square_double_gamma(d::Float64)

Calculate value of Buhnerkempe double gamma transmission kernel.\n
Accepts distance^2 instead of distance.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function USDOS2_kernel_square_double_gamma(d::Float64)

    #List of required parameter values
    k1::Float64 = 0.11489682
    k2::Float64 = 1600.0
    k3::Float64 = 9.2

    #Express kernel formulation
    KernelVal::Float64 = k1 / (1 + ((d^(0.5*k3))/(k2^k3)))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end


#-------------------------------------------------------------------------------
### 4: Gaussian kernel
#-------------------------------------------------------------------------------

"""
    gaussian_kernel(d::Float64)

Bivariate normal transmission kernel.

Inputs:
- `d::Float64`: Squared euclidean distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function gaussian_kernel(d::Float64)

    #List of required parameter values
    L_1_norm::Float64 = 0.2*(1000.0*1000.0)
    length_scale::Float64 = 3.0*1000.0
    fact_1::Float64 = 1.0 / (2*pi*length_scale*length_scale)

    #Express kernel formulation
    KernelVal::Float64 = L_1_norm * fact_1 * exp(-(d / (2*length_scale*length_scale)))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end

#-------------------------------------------------------------------------------
### 5: Brand a3 kernel
#-------------------------------------------------------------------------------

"""
    Construct_brand_kernel_a3(MaxDist::Int64)

Create lookup vector for Brand a3 transmission kernel. Entry for each one metre increment.

Inputs:
- `MaxDist::Int64`: Maximum distance to compute value of transmission kernel for.

Outputs:
- `KernelLookUpVec::Array{Float64,1}`: Kernel term, with entry for each one metre increment.

Location: SpatialKernels.jl
"""
function Construct_brand_kernel_a3(MaxDist::Int64)

    #Create lookup vector. Entry for each one metre increment

    #Initialise lookup array
    KernelLookUpVec = zeros(Float64,MaxDist)

    #Iterate over desired distances. Assign kernel value to array
    for DistItr = 1:MaxDist
        KernelLookUpVec[DistItr] = brand_kernel_a3(DistItr)
    end

    return KernelLookUpVec::Array{Float64,1}
end

"""
    brand_kernel_a3(d::Float64)

Calculate value of Brand a3 transmission kernel.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function brand_kernel_a3(d::Float64)

    #Transform distance to km
    d = d*0.001

    #Note, we assume distance is input in metres format!

    #List of required parameter values
    beta::Float64 = 0.115
    L3::Float64 = 10.0
    N3::Float64 = 1.0 / (2*pi)

    #Express kernel formulation
    KernelVal::Float64 = beta * N3 * L3 * ((L3^2 + d^2)^-1.5)

    #Return kernal value for distance d from function
    return KernelVal::Float64
end

#-------------------------------------------------------------------------------
### 6: Brand a4 kernel
#-------------------------------------------------------------------------------
"""
    brand_kernel_a4(d::Float64)

Calculate value of Brand a4 transmission kernel.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function brand_kernel_a4(d::Float64)

    #Transform distance to km
    #d = d*0.001

    #Note, we assume distance is input in km format!

    #List of required parameter values
    beta::Float64 = 0.115
    L4::Float64 = 20.0 / (0.5*pi)
    N4::Float64 = L4/pi

    #Express kernel formulation
    KernelVal::Float64 = beta * N4 * L4 * ((L4^2 + d^2)^-2.0)

    #Return kernal value for distance d from function
    return KernelVal::Float64
end

#-------------------------------------------------------------------------------
### 7: Brand a5 kernel
#-------------------------------------------------------------------------------
"""
    brand_kernel_a5(d::Float64)

Calculate value of Brand a5 transmission kernel.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function brand_kernel_a5(d::Float64)

    #Transform distance to km
    #d = d*0.001

    #Note, we assume distance is input in km format!

    #List of required parameter values
    beta::Float64 = 0.115
    L5::Float64 = 20.0
    N5::Float64 = (3.0 * L5 * L5) / (2 * pi)

    #Express kernel formulation
    KernelVal::Float64 = beta * N5 * L5 * ((L5^2 + d^2)^-2.5)

    #Return kernal value for distance d from function
    return KernelVal::Float64
end

#-------------------------------------------------------------------------------
### 8: Hayama (from Hayama et al 2015)
#-------------------------------------------------------------------------------
"""
    hayama_kernel(d::Float64)

Calculate value of Hayama transmission kernel (from Hayama et al 2015).

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function hayama_kernel(d::Float64)

    #Transform distance to km
    #d = d / 1000.0

    #Note, we assume distance is input in km format!

    #List of required parameter values
    Hayama_r0::Float64 = 0.58
    Hayama_h0::Float64 = 0.00074
    Hayama_alpha::Float64 = 2.47

    #Express kernel formulation
    KernelVal::Float64 = Hayama_h0 * ( (1 + (d/Hayama_r0))^(-Hayama_alpha))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end


#-------------------------------------------------------------------------------
### 9: USDOSv2 (from Tsao et al 2020)
#-------------------------------------------------------------------------------

"""
    Construct_USDOS2_kernel_v2(MaxDist::Int64)

Create lookup vector for USDOSv2 transmission kernel (from Tsao et al 2020). Entry for each one metre increment.

Inputs:
- `MaxDist::Int64`: Maximum distance to compute value of transmission kernel for.

Outputs:
- `KernelLookUpVec::Array{Float64,1}`: Kernel term, with entry for each one metre increment.

Location: SpatialKernels.jl
"""
function Construct_USDOS2_kernel_v2(MaxDist::Int64)

    #Create lookup vector. Entry for each one metre increment

    #Initialise lookup array
    KernelLookUpVec = zeros(Float64,MaxDist)

    #Iterate over desired distances. Assign kernel value to array
    for DistItr = 1:MaxDist
        KernelLookUpVec[DistItr] = USDOS2_kernel_v2(DistItr)
    end

    return KernelLookUpVec::Array{Float64,1}
end

"""
    USDOS2_kernel_v2(d::Float64)

Calculate value of USDOSv2 transmission kernel (from Tsao et al 2020).

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function USDOS2_kernel_v2(d::Float64)

    #List of required parameter values
    USDOS_v2_k1::Float64 = 1.46e-08 #1.293833e-08
    USDOS_v2_k2::Float64 = 1686.16 #2116.798
    USDOS_v2_k3::Float64 = 2.267 #2.38

    #Express kernel formulation
    KernelVal::Float64 = USDOS_v2_k1 / (1 + ((d/USDOS_v2_k2)^USDOS_v2_k3))

    #Return kernal value for distance d from function
    return KernelVal::Float64
end

#-------------------------------------------------------------------------------
### 10: GB FMD kernel (from Keeling & Rohani)
#-------------------------------------------------------------------------------

"""
    Construct_GB_FMD_kernel(MaxDist::Int64)

Create lookup vector for FMD transmission kernel (from Keeling & Rohani 2008). Entry for each one metre increment.

Inputs:
- `MaxDist::Int64`: Maximum distance to compute value of transmission kernel for.

Outputs:
- `KernelLookUpVec::Array{Float64,1}`: Kernel term, with entry for each one metre increment.

Location: SpatialKernels.jl
"""
function Construct_GB_FMD_kernel(MaxDist::Int64)

    #Create lookup vector. Entry for each one metre increment

    #Initialise lookup array
    KernelLookUpVec = zeros(Float64,MaxDist)

    #Iterate over desired distances. Assign kernel value to array
    MaxItr = min(MaxDist,60000)
    for DistItr = 1:MaxItr
        if DistItr <= 117 #Distances 117m and below
            KernelLookUpVec[DistItr] = 0.3093
        elseif DistItr <= 60000  #Distance up to 60km
            KernelLookUpVec[DistItr] = GB_FMD_kernel(DistItr)
        end
        #Note, all distance above 60km are zero.
        #Kernel values at these distances were already set to zero when array was initialised.
    end


    return KernelLookUpVec::Array{Float64,1}
end

"""
    GB_FMD_kernel(d::Float64)

Calculate value of FMD transmission kernel (from Keeling & Rohani 2008).

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `KernelVal::Float64`: Value of transmission kernel at distance d.

Location: SpatialKernels.jl
"""
function GB_FMD_kernel(d::Float64)

    #Convert distance to km
    d_km = d*0.001

    #Get squared distance
    d_squared = d_km^2

    #Initialise KernelVal
    KernelVal = 0.

    #Compute kernel value. Enforce threshold values for distances under 0.1km and above 60km
    if d_squared < 0.0138
        KernelVal = 0.3093
    elseif d_squared < 3600. #distance between 0.1km to 60km

        #@evalpoly(z, c...)
        #Evaluate the polynomial ∑(c[k]*z^{k−1}) for the coefficients c[1], c[2], ...;
        #that is, the coefficients are given in ascending order by power of z.
        #This macro expands to efficient inline code that uses either Horner's method or,
        #for complex z, a more efficient Goertzel-like algorithm.

        KernelVal = exp(@evalpoly(d_squared, -3.231772, -0.609262, -1.30519e-1, -3.3687e-2, 3.3966e-3, 9.5628e-4, -9.2123e-5))

        # @time KernelVal = exp( (-9.2123e-5*d_squared*d_squared*d_squared*d_squared*d_squared*d_squared) +
        #                                         (9.5628e-4*d_squared*d_squared*d_squared*d_squared*d_squared) +
        #                                           (3.3966e-3*d_squared*d_squared*d_squared*d_squared) -
        #                                          (3.3687e-2*d_squared*d_squared*d_squared) -
        #                                          (1.30519e-1*d_squared^2) -
        #                                          (0.609262*d_squared) - 3.231772)
    end

    #Return kernal value for distance d from function
    return KernelVal::Float64
end
