
# DGP Julia project (JuliaDGP.jl)
#
# some utilities
#
# AM

# collecting the common elements from a DGP and Realization instances
function common_elements(dgp::DGP{T},R::Realization{T}) where {T}
   if !are_compatible(dgp,R) throw(ArgumentError("The two instances differ in dimensionality")) end

   # interating and collecting
   el = Vector{T}();
   for v in keys(R.coords)
      if v in dgp.order
         if !(v in el)
            push!(el,v);
         end
      end
   end

   # the vector is ready
   return el;
end

# as above... the order of the arguments is not important
function common_elements(R::Realization{T},dgp::DGP{T}) where {T}
   return common_elements(dgp,R);
end

# counting the number of common elements from a DGP and Realization instances
function count_common_elements(dgp::DGP{T},R::Realization{T}) where {T}
   return length(common_elements(dgp,R));
end

# as above... the order of the arguments is not important
function count_common_elements(R::Realization{T},dgp::DGP{T}) where {T}
   return length(common_elements(R,dgp));
end

# extracting the amino acid sequence from a Vector{Atom} instance
function sequence(atoms::Vector{Atom})
   if (length(atoms) == 0) throw(ArgumentError("Input Vector{Atom} has zero length")) end

   # extracting the sequence
   seq = Dict{Int64,Char}();
   for atom in atoms
      if (!haskey(seq,atom.resnum)) push!(seq,atom.resnum => aa_name_convert(atom.resname)) end
   end

   # reordering by residue number
   resnums = collect(keys(seq));
   sort!(resnums);
   primary = "";
   for k in resnums
      primary = primary * seq[k];
   end

   # the String with the amino acid sequence is ready
   return primary;
end

# extracting the amino acid sequence from a PDB file
function sequence(PDBfile::String)
   ex = get_extension(PDBfile);
   if ex != "pdb"
      throw(ArgumentError("Input file is supposed to be a PDB file; check if extension is coherent with content"));
   end

   # loading the PDB file
   atoms = PDBTools.readPDB(PDBfile);
   if (atoms == nothing || length(atoms) == 0)  throw(ArgumentError("Something went wrong when reading the PDB file")) end

   # extracting the sequence and returning
   return sequence(atoms);
end

# dot conversion for the undirected graphs (SimpleGraphs)
function dot(g::UndirectedGraph{T}) where {T}
   dot_string = "graph {\n";
   for (u,v) in g.E
       dot_string = dot_string * "  ";
       dot_string = dot_string * "\"" * string(u) * "\" -- \"" * string(v) * "\";\n"
   end
   return dot_string * "}";
end

# extracting the extension of a file
function get_extension(filename::String)
   fileparts = split(filename,".");
   nparts = length(fileparts);
   extension = String(fileparts[nparts]);
   return extension;
end

