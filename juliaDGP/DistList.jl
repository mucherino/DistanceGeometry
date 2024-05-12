
# DGP Julia project (JuliaDGP.jl)
#
# DistList abstract type
#
# every structure extending this abstract type is supposed
# to have at least these two attributes:
# - K::Int64 for the dimensionality
# - distlist::Dict{Tuple{T,T},Distance} for the list of distances
#
# includes Dimensionality and Distance
# using SimpleGraphs
#
# AM

abstract type DistList{T} <: Dimensionality end

# counting the number of distances
function nd(dl::DistList{T}) where {T}
   return length(keys(dl.distlist));
end

# counting the number of exact distances
function nexact(dl::DistList{T}) where {T}
   n = 0;
   for pair in keys(dl.distlist)
      D = dl.distlist[pair];
      if is_exact(D)
         n = n + 1;
      end
   end
   return n;
end

# counting the number of interval distances
function ninterval(dl::DistList{T}) where {T}
   n = 0;
   for pair in keys(dl.distlist);
      D = dl.distlist[pair];
      if is_interval(D)
         n = n + 1;
      end
   end
   return n;
end

# extracting the elements from a DistList type
function elements(dl::DistList{T}) where {T}
   el = Vector{T}();
   for (u,v) in keys(dl.distlist)
      if !(u in el) push!(el,u) end
      if !(v in el) push!(el,v) end
   end
   return el;
end

# counting the number of elements in the DistList type
function nelements(dl::DistList{T}) where {T}
   return length(elements(dl));
end

# showing the details of the DistList instance
function details(dl::DistList{T}) where {T}
   for (u,v) in keys(dl.distlist)
      println("(",u,",",v,") => ",value(dl.distlist[(u,v)]));
   end
end

# collecting the distances from a given DistList instance related to a subset of its elements
function sublist(dl::DistList{T},subset::Vector{T}) where {T}
   el = elements(dl);
   for v in subset
      if !(v in el) throw(ArgumentError("Elements of the given subset are not in the DistList instance")) end
   end

   # collecting the distances
   distvect = Vector{Distance}();
   for u in subset
      for v in subset
         if u != v
            if haskey(dl.distlist,(u,v))
               push!(distvect,dl.distlist[(u,v)]);
            end
         end
      end
   end

   # the vector of distances is ready
   return distvect;
end

# verifying whether a given set of T elements form a clique in DistList instance
function is_clique(dl::DistList{T},clique::Vector{T}) where {T}
   n = length(clique);
   if (n == 0) return true end

   # verifying whether the given T elements belong to the DistList instance
   el = elements(dl);
   for k = 1:n
      if !(clique[k] in el) return false end
   end
   if (n == 1) return true end  # trivial clique

   # verifying whether the elements form a nontrivial clique
   k = 2;
   for i = 1:n
      u = clique[i];
      for j in 1:i - 1
         v = clique[j];
         if !haskey(dl.distlist,(u,v)) && !haskey(dl.distlist,(v,u)) return false end
      end
   end

   # at this point we have a clique
   return true;
end

# defining the undirected graph related to the DistList structure
function graph(dl::DistList{T}) where {T}
   g = UndirectedGraph{T}();

   # we construct the edge set
   for (u,v) in keys(dl.distlist)
      add!(g,u,v);
   end

   # the graph is ready
   return g;
end

