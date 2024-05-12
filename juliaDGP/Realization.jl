
# DGP Julia project (see JuliaDGP.jl)
#
# Realization structure
#
# importing Base.show
# using DataStructures
# including Dimensionality
#
# AM

struct Realization{T} <: Dimensionality
   K::Int64;  # dimension
   N::Int64;  # length
   coords::OrderedDict{T,Vector{Float64}};  # map of T types to coordinate vectors

   # constructor for random Realization instance
   function Realization(K::Int64,N::Int64)
      if (K <= 0) throw(ArgumentError("The dimension cannot be nonpositive")) end
      if (N <= 0) throw(ArgumentError("The size of the DGP instance cannot be nonpositive")) end

      # generating random coordinate vectors and associating them to the first N integer numbers
      coords = OrderedDict{Int64,Vector{Float64}}();
      for i in 1:N
         push!(coords,i => rand(Float64,K));
      end

      # the Realization is ready
      new{Int64}(K,N,coords);
   end

   # constructor from a list of Atom instances (PDBTools.Atom)
   function Realization(atoms::Vector{Atom})
      if (length(atoms) == 0) throw(ArgumentError("Input Vector{Atom} has zero length")) end

      # constructing the Realization instance
      coords = OrderedDict{myAtom,Vector{Float64}}();
      for atom in atoms
         myatom = myAtom(atom);
         push!(coords,myatom => [atom.x,atom.y,atom.z]);
      end

      # the Realization is ready
      new{myAtom}(3,length(atoms),coords);
   end

   # constructor from PDB file
   function Realization(PDBfile::String,modelId::Int64)
      ex = get_extension(PDBfile);
      if ex != "pdb"
         throw(ArgumentError("Input file is supposed to be a PDB file; check if extension is coherent with content"));
      end
      if (modelId <= 0) throw(ArgumentError("Model ID cannot be nonpositive")) end

      # loading the PDB file
      atoms = PDBTools.readPDB(PDBfile);
      if (atoms == nothing || length(atoms) == 0)  throw(ArgumentError("Something went wrong while reading the PDB file")) end
      atoms = PDBTools.select(atoms,"model = " * string(modelId));
      if (atoms == nothing || length(atoms) == 0)  throw(ArgumentError("Something went wrong while selecting the requested model")) end

      # invoking the constructor from the list of Atom instances
      Realization(atoms);
   end

   # constructor from PDB file (always first model)
   function Realization(PDBfile::String)
      Realization(PDBfile,1);
   end

   # overriding Base show function
   function Base.show(io::IO,R::Realization{T}) where {T}
      print(io,"Realization{",nameof(T),"} (",R.N," elements in dimension ",R.K,")");
   end
end

# length of a Realization instance
function length(R::Realization{T}) where {T}
   return R.N;
end

# extracting the elements from a Realization instance
function elements(R::Realization{T}) where {T}
   return collect(keys(R.coords));
end

# showing the details of a Realization instance
function details(R::Realization{T}) where {T}
   println("Realization{",nameof(T),"} (K = ",R.K,") {");
   for v in keys(R.coords)
      println("  ",v," => ",R.coords[v]);
   end
   print("}");
end

