
# DGP Julia project (see JuliaDGP.jl)
#
# ErrorList type
#
# importing Base.show
# including DistList, Realization, DGP
#
# AM

struct ErrorList{T} <: DistList{T}
   K::Int64;  # dimension
   distlist::Dict{Tuple{T,T},Distance};  # list of errors on distances

   # constructor from a DGP instance and one possible Realization instance (a possible solution)
   function ErrorList(dgp::DGP{T},R::Realization{T}) where {T}
      distlist = Dict{Tuple{T,T},Distance}();

      # computing and verifying the distances
      for (u,v) in keys(dgp.distlist)
         if haskey(R.coords,u) && haskey(R.coords,v)
            xu = R.coords[u];
            xv = R.coords[v];
            computed = Distance(xu,xv);
            expected = dgp.distlist[(u,v)];
            error = distance(computed,expected);
            if error > 0.0
               push!(distlist,(u,v) => Distance(error));
            end
         end
      end

      # the error list is ready
      new{T}(dgp.K,distlist);
   end

   # the same constructor as above... the order of the arguments is not important
   function ErrorList(R::Realization{T},dgp::DGP{T}) where {T}
      ErrorList(dgp,R);
   end

   # overriding Base show function
   function Base.show(io::IO,e::ErrorList{T}) where {T}
      n = nelements(e);
      m = nd(e);
      print(io,"ErrorList{",nameof(T),"} (K = ",e.K," with ",n," elements and ",m," errors)");
   end
end

# looking for the smallest error
function minerror(e::ErrorList{T}) where {T}
   smallest = nothing;
   for (u,v) in keys(e.distlist)
      current = value(e.distlist[(u,v)]);
      if smallest == nothing 
         smallest = ((u,v) => current);
      elseif current < smallest.second
         smallest = ((u,v) => current);
      end
   end
   return smallest;
end

# looking for the largest error
function maxerror(e::ErrorList{T}) where {T}
   largest = nothing;
   for (u,v) in keys(e.distlist)
      current = value(e.distlist[(u,v)]);
      if largest == nothing
         largest = ((u,v) => current);
      elseif largest.second < current
         largest = ((u,v) => current);
      end
   end
   return largest;
end

# computing the mde
function mde(e::ErrorList{T}) where {T}
   error = 0.0;
   for (u,v) in keys(e.distlist)
      error = error + value(e.distlist[(u,v)]);
   end
   return error / nd(e);
end

