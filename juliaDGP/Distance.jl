
# DGP Julia project (JuliaDGP.jl)
#
# Distance type
#
# importing Base.show and PDBTools.distance
#
# AM

struct Distance
   lb::Float64;
   ub::Float64;
   uncertainty::UInt8;

   # generic constructor
   function Distance(l::Float64,u::Float64,un::Int64)
      if (l > u) throw(ArgumentError("Lower bound is larger than upper bound")) end
      if (un < 0) throw(ArgumentError("The uncertainty level cannot be negative")) end
      new(l,u,un);
   end

   # generic constructor with Bool type
   function Distance(l::Float64,u::Float64,un::Bool)
      if (l > u) throw(ArgumentError("Lower bound is larger than upper bound")) end
      uncertainty = 0;
      if (un) uncertainty = 1 end
      new(l,u,uncertainty);
   end

   # constructor without uncertainty
   function Distance(l::Float64,u::Float64)
      if (l > u) throw(ArgumentError("Lower bound is larger than upper bound")) end
      new(l,u,0);
   end

   # constructor for exact distances
   function Distance(d::Float64)
      new(d,d,0);
   end

   # constructor from a pair of Vector types
   function Distance(u::Vector{Float64},v::Vector{Float64})
      K = length(u)
      if (K != length(v)) throw(ArgumentError("The two Vector instances need to have the same length")) end

      # computing the Euclidean distance
      value = 0.0;
      for i = 1:K
         diff = u[i] - v[i];
         diff = diff*diff;
         value = value + diff;
      end
      value = sqrt(value);

      # new distance
      new(value,value,0);
   end

   # overriding equality operator '=='
   function Base.:(==)(d1::Distance,d2::Distance)
      if d1.lb != d2.lb
         return false;
      elseif d1.ub != d2.ub
         return false;
      elseif d1.uncertainty == 0
         if d2.uncertainty != 0 return false end
      else
         if d2.uncertainty == 0 return false end
      end
      return true;
   end

   # overriding hash function
   function Base.hash(d::Distance)
      h = hash(d.lb) + hash(d.ub);
      if (d.uncertainty > 0) h = h << 1 end
      return h;
   end

   # overriding Base show function
   function Base.show(io::IO,d::Distance)
      print(io,"Distance(",d.lb)
      if d.lb != d.ub
         print(io,",",d.ub);
      end
      if d.uncertainty == 0
         print(io,")");
      else
         print(io,",uncertainty level ",d.uncertainty,")");
      end
   end
end

# is the Distance instance holding an exact distance?
function is_exact(d::Distance)
   return d.lb == d.ub;
end

# are the Distance instances in the Vector all holding exact distances?
function are_exact(D::Vector{Distance})
   for d in D
      if !is_exact(d) return false end
   end
   return true;
end

# counting the number of Distance instances in the Vector that are exact
function nexact(D::Vector{Distance})
   count = 0;
   for d in D
      if (is_exact(d)) count = count + 1 end
   end
   return count;
end

# is the Distance instance holding an interval distance?
function is_interval(d::Distance)
   return d.lb != d.ub;
end

# are the Distance instances in the Vector all holding exact distances?
function are_intervals(D::Vector{Distance})
   for d in D
      if !is_interval(d) return false end
   end
   return true;
end

# counting the number of Distance instances in the Vector that are intervals
function ninterval(D::Vector{Distance})
   count = 0;
   for d in D
      if (is_interval(d)) count = count + 1 end
   end
   return count;
end

# is the Distance instance uncertain?
function is_uncertain(d::Distance)
   return d.uncertainty > 0;
end

# are the Distance instances in the Vector all holding uncertain distances?
function are_uncertain(D::Vector{Distance})
   for d in D
      if !is_uncertain(d) return false end
   end
   return true;
end

# counting the number of Distance instances in the Vector that are uncertain
function nuncertain(D::Vector{Distance})
   count = 0;
   for d in D
      if (is_uncertain(d)) count = count + 1 end
   end
   return count;
end

# computing the range of the Distance instance
function range(d::Distance)
   return d.ub - d.lb;
end

# computing the (mean) value of the Distance instance
function value(d::Distance)
   if is_exact(d) return d.lb end
   v = 0.5*(d.lb + d.ub);
   return v;
end

# generating the intersection of two Distance instances
function intersect(d1::Distance,d2::Distance)
   if is_uncertain(d1) || is_uncertain(d2)
      throw(ArgumentError("It is not possible to intersect uncertain distances"));
   end
   t1 = is_exact(d1);
   t2 = is_exact(d2);
   if t1 && t2  # both are exact
      if d1.lb == d2.lb
         return Distance(d1.lb,d1.lb);
      end
   elseif t1 && !t2  # the second is an interval
      if d2.lb <= d1.lb && d1.lb <= d2.ub
         return Distance(d1.lb,d1.lb);
      end
   elseif !t1 && t2  # the first is an interval
      if d1.lb <= d2.lb && d2.lb <= d1.ub
         return Distance(d2.lb,d2.lb)
      end
   else  # both are intervals
      i1 = d1;
      i2 = d2;
      if i1.lb > i2.lb
         i1 = d2;
         i2 = d1;
      end
      if i1.ub >= i2.lb
         return Distance(i2.lb,i1.ub);
      end
   end
   return nothing; 
end

# computing the "distance" between two Distance instances
function distance(d1::Distance,d2::Distance)
   if intersect(d1,d2) != nothing
      return 0.0;
   end
   return min(abs(d1.lb - d2.ub),abs(d2.lb - d1.ub));
end

# creating a new Distance from a previous one but with increased uncertainty
# (the hash code is the same of the original instance)
function with_higher_uncertainty(d::Distance)
   if (d.uncertainty == 0) throw(ArgumentError("This Distance type does not hold an uncertain distance")) end
   return Distance(d.lb,d.ub,d.uncertainty+1);
end

