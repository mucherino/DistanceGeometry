
# DGP Julia project (see JuliaDGP.jl)
#
# myAtom type
#
# importing Base.show and using PDBTools
#
# AM

struct myAtom
   name::String;
   resnum::Int64;
   resname::String;

   # basic constructor with some argument verifications
   function myAtom(name::String,resnum::Int64,resname::String)
      # name
      if (length(name) == 0) throw(ArgumentError("myAtom name cannot be empty")) end
      if (!isletter(name[1])) throw(ArgumentError("myAtom name should always start with a letter")) end
      i = 1;
      while (i <= length(name) && isletter(name[i])) i = i + 1 end
      while i <= length(name)
         invalid = ArgumentError("Valid myAtom names should begin with letters, then digits, possibly followed by one ending wildcard #");
         if (i < length(name) && !isnumeric(name[i])) throw(invalid) end
         if (i == length(name) && !isnumeric(name[i]) && name[i] != '#') throw(invalid) end
         i = i + 1;
      end

      # resnum
      if (resnum <= 0) throw(ArgumentError("myAtom resnum cannot be nonpositive")) end

      # resname
      converted = aa_name_convert(resname);
      if (converted == nothing) throw(ArgumentError("The given residue name is not a valid name")) end
      if (length(resname) != 3) resname = converted end
     
      # myAtom instance is ready
      new(name,resnum,resname);
   end

   # constructor from an Atom instance from PDBTools
   function myAtom(atom::Atom)
      new(atom.name,atom.resnum,atom.resname);
   end

   # overriding equality operator
   function Base.:(==)(a::myAtom,b::myAtom)
      return a.name == b.name && a.resnum == b.resnum && a.resname == b.resname;
   end

   # overriding hash code function
   function Base.hash(a::myAtom)
      return hash(a.name) + hash(a.resnum) + hash(a.resname);
   end

   # overriding string conversion function
   function Base.string(a::myAtom)
      return a.resname * string(a.resnum) * "-" * a.name;
   end

   # overriding show function
   function Base.show(io::IO,a::myAtom)
      print(io,string(a));
   end
end

# verifying whether two myAtom instances are equivalent
function are_equivalent(a::myAtom,b::myAtom)
   if (a.resnum != b.resnum) return false end
   if (a.resname != b.resname) return false end
   alen = length(a.name);
   blen = length(b.name);
   if (alen != blen) return false end
   if a.name[alen] == '#' || b.name[blen] == '#'
      if (a.name[1:alen-1] != b.name[1:blen-1]) return false end
   end
   return true;
end

# converting from Vector{Atom} to Vector{myAtom}
function myAtoms(atoms::Vector{Atom})
   myatoms = Vector{myAtom}();
   for a in atoms
      push!(myatoms,myAtom(a));
   end
   return myatoms;
end

# converting from myAtom instance to PDBTools.Atom instance
function PDBToolAtom(index::Int64,a::myAtom)
   if (index < 0) throw(ArgumentError("The PDBTools.Atom index cannot be negative")) end
   type = a.name[1:1];
   return PDBTools.Atom(index,index,a.name,a.resname,"A",a.resnum,a.resnum,0.0,0.0,0.0,0.0,0.0,1,"-",type,nothing,Dict{Symbol,Any}());
end

