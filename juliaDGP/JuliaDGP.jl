
# DGP Julia project
#
# AM

import Base.show
import Base.length
import Base.string
import PDBTools.distance
using PDBTools
using SimpleGraphs
using Bijections
using DataStructures
using Printf
include("Distance.jl")
include("myAtom.jl")
include("Aminoacids.jl")
include("Dimensionality.jl")
include("DistList.jl")
include("Realization.jl")
include("DGP.jl")
include("ErrorList.jl")
include("Orders.jl")
include("Utils.jl")
name = PDBTools.name

