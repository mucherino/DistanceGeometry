
# DGP Julia project (see JuliaDGP.jl)
#
# DGP type
#
# importing Base.show
# using DataStructures
# including DistList, Orders, Realization, myAtom, Utils
#
# AM

struct DGP{T} <: DistList{T}
   K::Int64;  # dimension
   distlist::Dict{Tuple{T,T},Distance};  # list of distances
   order::Vector{T};  # order for the elements of type T

   # constructor for random fully-connected instance from given Realization
   function DGP(R::Realization{T}) where {T}
      order = vertex_order(R);

      # computing all the distances
      n = length(order);
      distlist = Dict{Tuple{T,T},Distance}();
      for i in 1:n
         u = order[i];
         for j in 1:i - 1
            v = order[j];
            push!(distlist,(u,v) => Distance(R.coords[u],R.coords[v]));
         end
      end

      # DGP instance is ready
      new{T}(R.K,distlist,order);
   end

   # constructor from a given Realization instance and distance cut-off
   function DGP(R::Realization{T},threshold::Float64) where {T}
      if (threshold <= 0.0) throw(ArgumentError("The specified threshold is 0.0 or even negative")) end
      order = vertex_order(R);

      # computing and selecting distances
      n = length(order);
      distlist = Dict{Tuple{T,T},Distance}();
      for i = 1:n
         u = order[i];
         for j in 1:i - 1
            v = order[j];
            dist = Distance(R.coords[u],R.coords[v]);
	    if dist.ub <= threshold
               push!(distlist,(u,v) => dist);
            end
         end
      end

      # DGP instance is ready
      new{T}(R.K,distlist,order);
   end

   # constructor for a paradoxical DGP from a given Realization instance
   # -> "cut" is the (constant) number of edges crossing DGP elements
   function DGP(R::Realization{T},cut::Int64) where {T}
      n = R.N;
      if (cut <= 0) throw(ArgumentError("THe specified cut value is 0 or even negative")) end
      if (cut > n) throw(ArgumentError("The value of cut cannot be larger than the Realization size")) end
      if (cut == n) return DGP(R) end  # it's going to be a complete DGP instance

      # preparing the vectors of elements and references
      elements = vertex_order(R);
      refs = Deque{T}();
      for i in n-cut+1:n
         push!(refs,elements[i]);
      end

      # constructing the paradoxical instance
      distlist = Dict{Tuple{T,T},Distance}();
      for i in 1:n
         v = elements[i];
         for u in refs
            push!(distlist,(u,v) => Distance(R.coords[u],R.coords[v]));
         end
         popfirst!(refs);
         push!(refs,v);
      end

      # paradoxical DGP is ready
      new{T}(R.K,distlist,elements);
   end

   # constructor from a given Realization instance and its connection graph
   function DGP(R::Realization{T},g::UndirectedGraph{T}) where {T}
      elements = vertex_order(R);  # initial order with *all* vertices

      # selecting only distances for which an edge is in the graph g
      n = length(elements);
      order = OrderedDict{Int64,T}();
      distlist = Dict{Tuple{T,T},Distance}();
      for i = 1:n
         u = elements[i];
         for j in 1:i - 1
            v = elements[j];
            if (u,v) in g.E || (v,u) in g.E
               push!(order,i => u);
               push!(order,j => v);
               push!(distlist,(u,v) => Distance(R.coords[u],R.coords[v]));
            end
         end
      end

      # DGP instance is ready
      new{T}(R.K,distlist,collect(values(order)));
   end

   # constructor from STAR file (NMR data, we ignore uncertain distances)
   function DGP(STARfile::String)
      ex = get_extension(STARfile);
      if ex != "str"
         throw(ArgumentError("Input file is supposed to be a STAR file; check if extension is coherent with content"));
      end

      # preparing the data structures
      order = Vector{myAtom}();
      distlist = Dict{Tuple{myAtom,myAtom},Distance}();

      # opening the STAR file
      star = open(STARfile);

      # reading the lines of the STAR file
      extracting = false;
      identifiers = Vector{String}();
      for line in readlines(star)
         if !extracting
            # looking for the blocks of data
            if contains(line,"loop_")
               extracting = true;
            end
         else
            if contains(line,"_Gen_dist_constraint.")
               # reading the identifiers
               splitline = split(line,'.');
               id = String(splitline[2]);
               push!(identifiers,id);
            elseif !isempty(identifiers)
               # reading the data if the searched identifiers were found
               n = length(identifiers);
               splitline = split(line);
               if n == length(splitline)
                  # preparing some local variables
                  resnum1 = nothing;  resnum2 = nothing;
                  resname1 = nothing; resname2 = nothing;
                  name1 = nothing;    name2 = nothing;
                  lb = nothing;       ub = nothing;
                  uncertain = false;

                  # reading
                  for i in 1:n
                     if contains(identifiers[i],"Member_logic_code")
                        if (String(splitline[i]) == "OR") uncertain = true; break end  # we skip it!
                     elseif contains(identifiers[i],"PDB_residue_no_1")
                        resnum1 = parse(Int64,String(splitline[i]));
                     elseif contains(identifiers[i],"PDB_residue_no_2")
                        resnum2 = parse(Int64,String(splitline[i]));
                     elseif contains(identifiers[i],"PDB_residue_name_1")
                        resname1 = String(splitline[i]);
                     elseif contains(identifiers[i],"PDB_residue_name_2")
                        resname2 = String(splitline[i]);
                     elseif contains(identifiers[i],"Auth_atom_ID_1")
                        name1 = String(splitline[i]);
                     elseif contains(identifiers[i],"Auth_atom_ID_2")
                        name2 = String(splitline[i]);
                     elseif contains(identifiers[i],"Distance_lower_bound_val")
                        lb = parse(Float64,String(splitline[i]));
                     elseif contains(identifiers[i],"Distance_upper_bound_val")
                        ub = parse(Float64,String(splitline[i]));
                     end
                  end

                  # defining the data structures
                  if !uncertain
                     if (resnum1 == nothing || resnum2 == nothing) throw(ArgumentError("Something went wrong...")) end
                     if (resname1 == nothing || resname2 == nothing) throw(ArgumentError("Something went wrong...")) end
                     if (name1 == nothing || name2 == nothing) throw(ArgumentError("Something went wrong...")) end
                     if (name1 == "HN")  name1 = "H" end
                     a = myAtom(name1,resnum1,resname1);
                     if !(a in order) push!(order,a) end
                     if (name2 == "HN")  name2 = "H" end
                     b = myAtom(name2,resnum2,resname2);
                     if !(b in order) push!(order,b) end
                     if (lb == nothing || ub == nothing) throw(ArgumentError("Something went wrong...")) end
                     D = Distance(lb,ub);
                     push!(distlist,(a,b) => D);
                  end
               end
            elseif contains(line,"stop_")
               # stopping
               identifiers = Vector{String}();
               extracting = false;
            end
         end
      end

      # closing the STAR file
      close(star);

      # DGP instance is ready
      new{myAtom}(3,distlist,order);
   end

   # overriding Base show function
   function Base.show(io::IO,dgp::DGP{T}) where {T}
      n = length(dgp.order);
      m = nd(dgp);
      print(io,"DGP{",nameof(T),"} (K = ",dgp.K," with ",n," elements and ",m," distances)");
   end
end

# showing the details of a DGP instance
function details(dgp::DGP{T}) where {T}
   println("DGP{",nameof(T),"} (K = ",dgp.K,") {");
   for u in dgp.order
      for v in dgp.order
         if u != v
            if haskey(dgp.distlist,(u,v))
               println("  (",u,",",v,") => ",dgp.distlist[(u,v)]);
            end
         end
      end
   end
   print("}");
end

# finding the references for a given element v of a DGP instance
function references(dgp::DGP{T},v::T) where {T}
   if !(v in dgp.order) throw(ArgumentError("The given element does not belong to the DGP instance")) end

   # looking for the position of the element in the order
   k = findfirst(x -> x == v,dgp.order);

   # looking for the references
   refs = Vector{T}();
   for i in 1:k - 1
      u = dgp.order[i];
      if haskey(dgp.distlist,(u,v)) || haskey(dgp.distlist,(v,u))
         push!(refs,u);
      end
   end

   # reference vector is ready
   return refs;
end

# printing in MDjeep format 
function MDjeep_format(dgp::DGP{T}) where {T}
   file = open("test_instance.nmr","w");
   len = length(dgp.order);
   for i in 1:len
      u = dgp.order[i];
      for j in 1:i - 1
         v = dgp.order[j];
         D = nothing;
         if (haskey(dgp.distlist,(u,v))) D = dgp.distlist[(u,v)] end
         if (haskey(dgp.distlist,(v,u))) D = dgp.distlist[(v,u)] end
         if D != nothing
            line = string(j) * " " * string(i) * " " * @sprintf("%20.16f",D.lb) * " " * @sprintf("%20.16f",D.ub);
	    if typeof(v) == myAtom && typeof(u) == myAtom
               line = line * " " * v.name * " " * u.name;
               line = line * " " * string(v.resnum) * " " * string(u.resnum);
               line = line * " " * v.resname * " " * u.resname;
            end
            line = line * "\n";
            write(file,line);
         end
      end
   end
   println("File 'test_instance.nmr' written in MDjeep format.");
   close(file);
end

