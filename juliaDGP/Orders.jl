
# DGP Julia project (JuliaDGP.jl)
#
# Order-related functions
#
# including DGP and Realization
#
# AM

# extracting the vertex order from a Realization instance
function vertex_order(R::Realization{T}) where {T}
   return collect(keys(R.coords));
end

# extracting the vertex order from a DGP instance
function vertex_order(dgp::DGP{T}) where {T}
   return copy(dgp.order);
end

