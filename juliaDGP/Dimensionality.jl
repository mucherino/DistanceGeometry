
# DGP Julia project (JuliaDGP.jl)
#
# Dimensionality abstract type
#
# every structure extending Dimensionality is supposed
# to have a K::Int64 attribute
#
# AM

abstract type Dimensionality end

# getting the dimension of the structure implementing this abstract type
function dimension(x::Dimensionality) :: Int64
   return x.K;
end

# verifying the compatibility in terms of dimensionality
function are_compatible(x::Dimensionality,y::Dimensionality) :: Bool
   return x.K == y.K;
end

