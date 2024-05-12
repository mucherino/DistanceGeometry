
# DGP Julia project (see JuliaDGP.jl)
#
# Tools for aminoacids
#
# using Bijection and SimpleGraphs
#
# AM

# bijection for amino acid codes
function codes_bijection()
   map = Bijection{Char,String}();
   push!(map,'A' => "ALA");
   push!(map,'R' => "ARG");
   push!(map,'N' => "ASN");
   push!(map,'D' => "ASP");
   push!(map,'B' => "ASX");  # ambiguous
   push!(map,'C' => "CYS");
   push!(map,'E' => "GLU");
   push!(map,'Q' => "GLN");
   push!(map,'Z' => "GLX");  # ambiguous
   push!(map,'G' => "GLY");
   push!(map,'H' => "HIS");
   push!(map,'I' => "ILE");
   push!(map,'L' => "LEU");
   push!(map,'K' => "LYS");
   push!(map,'M' => "MET");
   push!(map,'F' => "PHE");
   push!(map,'P' => "PRO");
   push!(map,'S' => "SER");
   push!(map,'T' => "THR");
   push!(map,'W' => "TRP");
   push!(map,'Y' => "TYR");
   push!(map,'V' => "VAL");
   return map;
end

# bijection for amino acid functions
function fun_dict()
   map = Dict{Char,String}();
   push!(map,'A' => "alanine");
   push!(map,'R' => "arginine");
   push!(map,'N' => "asparagine");
   push!(map,'D' => "aspartic");
   push!(map,'C' => "cysteine");
   push!(map,'E' => "glutamic");
   push!(map,'Q' => "glutamine");
   push!(map,'G' => "glycine");
   push!(map,'H' => "histidine");
   push!(map,'I' => "isoleucine");
   push!(map,'L' => "leucine");
   push!(map,'K' => "lysine");
   push!(map,'M' => "methionine");
   push!(map,'F' => "phenylalanine");
   push!(map,'P' => "proline");
   push!(map,'S' => "serine");
   push!(map,'T' => "threonine");
   push!(map,'W' => "tryptophan");
   push!(map,'Y' => "tyrosine");
   push!(map,'V' => "valine");
   return map;
end

# converting amino acid names
function aa_name_convert(name::String)
   map = codes_bijection();
   n = length(name);
   if n == 1
      # converting from 1-letter to 3-letter code
      c = name[1];
      if (haskey(map,c)) return map[c] end
      throw(ArgumentError("Unknown 1-letter code for amino acid"));
   elseif n == 3
      # converting from 3-letter to 1-letter code
      map = inv(map);
      if (haskey(map,name)) return map[name] end
      throw(ArgumentError("Unknown 3-letter code for amino acid"));
   else
      throw(ArgumentError("Amino acid codes can be composed either by 1 or 3 letters"));
   end
   return nothing;
end

# adding main backbone structure to existing graph
function backbone!(g::UndirectedGraph{myAtom},resnum::Int64,resname::String)
   N = myAtom("N",resnum,resname);
   H = myAtom("H",resnum,resname);
   CA = myAtom("CA",resnum,resname);
   C = myAtom("C",resnum,resname);
   O = myAtom("O",resnum,resname);

   # adding the bonds to the graph
   add!(g,N,H);
   add!(g,N,CA);
   add!(g,CA,C);
   add!(g,C,O);

   # we are done
   return nothing;
end

# main backbone graph only
function backbone(resnum::Int64,resname::String)
   g = UndirectedGraph{myAtom}();
   backbone!(g,resnum,resname);
   return g;
end

# adding alanine structure to existing graph
function alanine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "ALA";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB#",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,HB);

   # we are done
   return nothing;
end

# alanine graph only
function alanine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   alanine!(g,resnum);
   return g;
end

# adding arginine structure to existing graph
function arginine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "ARG";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   CD = myAtom("CD",resnum,resname);
   NE = myAtom("NE",resnum,resname);
   CZ = myAtom("CZ",resnum,resname);
   NH1 = myAtom("NH1",resnum,resname);
   NH2 = myAtom("NH2",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   HG = myAtom("HG#",resnum,resname);
   HD = myAtom("HD#",resnum,resname);
   HE = myAtom("HE",resnum,resname);
   HH1 = myAtom("HH1#",resnum,resname);
   HH2 = myAtom("HH2#",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD);
   add!(g,CD,NE);
   add!(g,NE,CZ);
   add!(g,CZ,NH1);
   add!(g,CZ,NH2);
   add!(g,CB,HB);
   add!(g,CG,HG);
   add!(g,CD,HD);
   add!(g,NE,HE);
   add!(g,NH1,HH1);
   add!(g,NH2,HH2);

   # we are done
   return nothing;
end

# arginine graph only
function arginine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   arginine!(g,resnum);
   return g;
end

# adding asparagine structure to existing graph
function asparagine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "ASN";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   OD1 = myAtom("OD1",resnum,resname);
   ND2 = myAtom("ND2",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   HD2 = myAtom("HD2#",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,OD1);
   add!(g,CG,ND2);
   add!(g,CB,HB);
   add!(g,ND2,HD2);

   # we are done
   return nothing;
end

# asparagine graph only
function asparagine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   asparagine!(g,resnum);
   return g;
end

# adding aspartic acid structure to existing graph
function aspartic!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "ASP";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   OD = myAtom("OD#",resnum,resname);
   HB = myAtom("HB",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,OD);
   add!(g,CB,HB);

   # we are done
   return nothing;
end

# aspartic acid graph only
function aspartic(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   aspartic!(g,resnum);
   return g;
end

# adding cysteine structure to existing graph
function cysteine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "CYS";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   SG = myAtom("SG",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   HG = myAtom("HG",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,SG);
   add!(g,CB,HB);
   add!(g,SG,HG);

   # we are done
   return nothing;
end

# cysteine graph only
function cysteine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   cysteine!(g,resnum);
   return g;
end

# adding glutamine structure to existing graph
function glutamine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "GLN";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   CD = myAtom("CD",resnum,resname);
   OE1 = myAtom("OE1",resnum,resname);
   NE2 = myAtom("NE2",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   HG = myAtom("HG#",resnum,resname);
   HE2 = myAtom("HE2",resnum,resname);
  
   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD);
   add!(g,CD,OE1);
   add!(g,CD,NE2);
   add!(g,CB,HB);
   add!(g,CG,HG);
   add!(g,NE2,HE2);

   # we are done
   return nothing;
end

# glutamine graph only
function glutamine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   glutamine!(g,resnum);
   return g;
end

# adding glutamic acid to existing graph
function glutamic!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "GLU";
   backbone!(g,resnum,resname);

   # prepating the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   CD = myAtom("CD",resnum,resname);
   OE = myAtom("OE",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   HG = myAtom("HG#",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD);
   add!(g,CD,OE);
   add!(g,CB,HB);
   add!(g,CG,HG);

   # we are done
   return nothing;
end

# glutamic acid graph only
function glutamic(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   glutamic!(g,resnum);
   return g;
end

# adding glycine structure to existing graph
function glycine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "GLY";
   backbone!(g,resnum,resname);
   add!(g,myAtom("CA",resnum,resname),myAtom("HA#",resnum,resname));
   return nothing;
end

# glycine graph only
function glycine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   glycine!(g,resnum);
   return g;
end

# adding histidine structure to existing graph
function histidine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "HIS";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   ND1 = myAtom("ND1",resnum,resname);
   CE1 = myAtom("CE1",resnum,resname);
   NE2 = myAtom("NE2",resnum,resname);
   CD2 = myAtom("CD2",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   HD1 = myAtom("HD1",resnum,resname);
   HE1 = myAtom("HE1",resnum,resname);
   HE2 = myAtom("HE2",resnum,resname);
   HD2 = myAtom("HD2",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,ND1);
   add!(g,ND1,CE1);
   add!(g,CE1,NE2);
   add!(g,NE2,CD2);
   add!(g,CD2,CG);
   add!(g,CB,HB);
   add!(g,ND1,HD1);
   add!(g,CE1,HE1);
   add!(g,NE2,HE2);
   add!(g,CD2,HD2);

   # we are done
   return nothing;
end

# histidine graph only
function histidine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   histidine!(g,resnum);
   return g;
end

# adding isoleucine structure to existing graph
function isoleucine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "ILE";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   CG1 = myAtom("CG1",resnum,resname);
   CG2 = myAtom("CG2",resnum,resname);
   CD1 = myAtom("CD1",resnum,resname);
   HG1 = myAtom("HG1#",resnum,resname);
   HG2 = myAtom("HG2#",resnum,resname);
   HD1 = myAtom("HD1#",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG1);
   add!(g,CB,CG2);
   add!(g,CG1,CD1);
   add!(g,CB,HB);
   add!(g,CG1,HG1);
   add!(g,CG2,HG2);
   add!(g,CD1,HD1);

   # we are done
   return nothing;
end

# isoleucine graph only
function isoleucine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   isoleucine!(g,resnum);
   return g;
end

# adding leucine structure to existing graph
function leucine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "LEU";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   HG = myAtom("HG",resnum,resname);
   CD1 = myAtom("CD1",resnum,resname);
   CD2 = myAtom("CD2",resnum,resname);
   HD1 = myAtom("HD1#",resnum,resname);
   HD2 = myAtom("HD2#",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD1);
   add!(g,CG,CD2);
   add!(g,CB,HB);
   add!(g,CG,HG);
   add!(g,CD1,HD1);
   add!(g,CD2,HD2);

   # we are done
   return nothing;
end

# leucine graph only
function leucine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   leucine!(g,resnum);
   return g;
end

# adding lysine structure to existing graph
function lysine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "LYS";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   HG = myAtom("HG#",resnum,resname);
   CD = myAtom("CD",resnum,resname);
   HD = myAtom("HD#",resnum,resname);
   CE = myAtom("CE",resnum,resname);
   HE = myAtom("HE#",resnum,resname);
   NZ = myAtom("NZ",resnum,resname);
   HZ = myAtom("HZ#",resnum,resname);

   # adding the new bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD);
   add!(g,CD,CE);
   add!(g,CE,NZ);
   add!(g,CB,HB);
   add!(g,CG,HG);
   add!(g,CD,HD);
   add!(g,CE,HE);
   add!(g,NZ,HZ);

   # we are done
   return nothing;
end

# lysine graph only
function lysine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   lysine!(g,resnum);
   return g;
end

# adding methionine structure to existing graph
function methionine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "MET";

   # preparing the backbone atoms
   N = myAtom("N",resnum,resname);
   H = myAtom("H#",resnum,resname);
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   C = myAtom("C",resnum,resname);
   O = myAtom("O",resnum,resname);

   # preparing the side chain atoms
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   SD = myAtom("SD",resnum,resname);
   CE = myAtom("CE",resnum,resname);
   HG = myAtom("HG#",resnum,resname);
   HE = myAtom("HE#",resnum,resname);

   # adding all bonds to the graph
   add!(g,N,H);
   add!(g,N,CA);
   add!(g,CA,C);
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,SD);
   add!(g,SD,CE);
   add!(g,CB,HB);
   add!(g,CG,HG);
   add!(g,CE,HE);
   add!(g,C,O);

   # we are done
   return nothing;
end

# methionine graph only
function methionine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   methionine!(g,resnum);
   return g;
end

# adding phenylalanine structure to existing graph
function phenylalanine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "PHE";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   CD1 = myAtom("CD1",resnum,resname);
   CD2 = myAtom("CD2",resnum,resname);
   CE1 = myAtom("CE1",resnum,resname);
   CE2 = myAtom("CE2",resnum,resname);
   CZ = myAtom("CZ",resnum,resname);
   HD1 = myAtom("HD1",resnum,resname);
   HD2 = myAtom("HD2",resnum,resname);
   HE1 = myAtom("HE1",resnum,resname);
   HE2 = myAtom("HE2",resnum,resname);
   HZ = myAtom("HZ",resnum,resname);

   # adding the bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD1);
   add!(g,CG,CD2);
   add!(g,CD1,CE1);
   add!(g,CD2,CE2);
   add!(g,CE1,CZ);
   add!(g,CE2,CZ);
   add!(g,CB,HB);
   add!(g,CD1,HD1);
   add!(g,CD2,HD2);
   add!(g,CE1,HE1);
   add!(g,CE2,HE2);
   add!(g,CZ,HZ);

   # we are done
   return nothing;
end

# phenylalanine graph only
function phenylalanine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   phenylalanine!(g,resnum);
   return g;
end

# adding proline structure to existing graph
function proline!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "PRO";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   N = myAtom("N",resnum,resname);
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB#",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   HG = myAtom("HG#",resnum,resname);
   CD = myAtom("CD",resnum,resname);
   HD = myAtom("HD",resnum,resname);

   # adding the bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD);
   add!(g,CD,N);
   add!(g,CB,HB);
   add!(g,CG,HG);
   add!(g,CD,HD);

   # we are done
   return nothing;
end

# proline graph only
function proline(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   proline!(g,resnum);
   return g;
end

# adding serine structure to existing graph
function serine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "SER";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB",resnum,resname);
   OG = myAtom("OG",resnum,resname);
   HG = myAtom("HG",resnum,resname);

   # adding the bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,OG);
   add!(g,CB,HB);
   add!(g,OG,HG);

   # we are done
   return nothing;
end

# serine graph only
function serine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   serine!(g,resnum);
   return g;
end

# adding threonine structure to existing graph
function threonine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "THR";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB",resnum,resname);
   OG1 = myAtom("OG1",resnum,resname);
   CG2 = myAtom("CG2",resnum,resname);
   HG1 = myAtom("HG1",resnum,resname);
   HG2 = myAtom("HG2#",resnum,resname);

   # adding the bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,OG1);
   add!(g,CB,CG2);
   add!(g,CB,HB);
   add!(g,OG1,HG1);
   add!(g,CG2,HG2);

   # we are done
   return nothing;
end

# serine graph only
function threonine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   threonine!(g,resnum);
   return g;
end

# adding tryptophan structure to existing graph
function tryptophan!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "TRP";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   CD1 = myAtom("CD1",resnum,resname);
   CD2 = myAtom("CD2",resnum,resname);
   NE1 = myAtom("NE1",resnum,resname);
   CE2 = myAtom("CE2",resnum,resname);
   CE3 = myAtom("CE3",resnum,resname);
   CZ2 = myAtom("CZ2",resnum,resname);
   CZ3 = myAtom("CZ3",resnum,resname);
   CH2 = myAtom("CH2",resnum,resname);
   HD1 = myAtom("HD1",resnum,resname);
   HE1 = myAtom("HE1",resnum,resname);
   HE3 = myAtom("HE3",resnum,resname);
   HZ2 = myAtom("HZ2",resnum,resname);
   HZ3 = myAtom("HZ3",resnum,resname);
   HH2 = myAtom("HH2",resnum,resname);

   # adding the bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD1);
   add!(g,CG,CD2);
   add!(g,CD1,NE1);
   add!(g,CD2,CE2);
   add!(g,NE1,CE2);
   add!(g,CD2,CE3);
   add!(g,CE2,CZ2);
   add!(g,CE3,CZ3);
   add!(g,CZ3,CH2);
   add!(g,CH2,CZ2);
   add!(g,CB,HB);
   add!(g,CD1,HD1);
   add!(g,NE1,HE1);
   add!(g,CE3,HE3);
   add!(g,CZ3,HZ3);
   add!(g,CH2,HH2);
   add!(g,CZ2,HZ2);

   # we are done
   return nothing;
end

# tryptophan graph only
function tryptophan(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   tryptophan!(g,resnum);
   return g;
end

# adding tyrosine structure to existing graph
function tyrosine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "TYR";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB",resnum,resname);
   CG = myAtom("CG",resnum,resname);
   CD1 = myAtom("CD1",resnum,resname);
   CD2 = myAtom("CD2",resnum,resname);
   CE1 = myAtom("CE1",resnum,resname);
   CE2 = myAtom("CE2",resnum,resname);
   CZ = myAtom("CZ",resnum,resname);
   OH = myAtom("OH",resnum,resname);
   HD1 = myAtom("HD1",resnum,resname);
   HD2 = myAtom("HD2",resnum,resname);
   HE1 = myAtom("HE1",resnum,resname);
   HE2 = myAtom("HE2",resnum,resname);
   HH = myAtom("HH",resnum,resname);

   # adding the bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG);
   add!(g,CG,CD1);
   add!(g,CG,CD2);
   add!(g,CD1,CE1);
   add!(g,CD2,CE2);
   add!(g,CE1,CZ);
   add!(g,CE2,CZ);
   add!(g,CZ,OH);
   add!(g,CB,HB);
   add!(g,CD1,HD1);
   add!(g,CD2,HD2);
   add!(g,CE1,HE1);
   add!(g,CE2,HE2);
   add!(g,OH,HH);

   # we are done
   return nothing;
end

# tyrosine graph only
function tyrosine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   tyrosine!(g,resnum);
   return g;
end

# adding valine structure to existing graph
function valine!(g::UndirectedGraph{myAtom},resnum::Int64)
   resname = "VAL";
   backbone!(g,resnum,resname);

   # preparing the necessary atoms
   CA = myAtom("CA",resnum,resname);
   HA = myAtom("HA",resnum,resname);
   CB = myAtom("CB",resnum,resname);
   HB = myAtom("HB",resnum,resname);
   CG1 = myAtom("CG1",resnum,resname);
   CG2 = myAtom("CG2",resnum,resname);
   HG1 = myAtom("HG1",resnum,resname);
   HG2 = myAtom("HG2",resnum,resname);

   # adding the bonds to the graph
   add!(g,CA,HA);
   add!(g,CA,CB);
   add!(g,CB,CG1);
   add!(g,CB,CG2);
   add!(g,CB,HB);
   add!(g,CG1,HG1);
   add!(g,CG2,HG2);

   # we are done
   return nothing;
end

# tyrosine graph only
function valine(resnum::Int64)
   g = UndirectedGraph{myAtom}();
   valine!(g,resnum);
   return g;
end

# constructing the protein graph from a given amino acid sequence
function protein(sequence::String)
   n = length(sequence)
   if (n == 0) throw(ArgumentError("The protein sequence has zero length")) end

   # loading the mapping from residue code to functions
   map = fun_dict();

   # constructing the first residue
   rescode = sequence[1];
   g = getfield(Main,Symbol(map[rescode]))(1);

   # constructing the remaining protein graph
   for i in 2:n
      prev = rescode;
      rescode = sequence[i];

      # adding the ith residue
      getfield(Main,Symbol(map[rescode]*"!"))(g,i);

      # including the peptide bond
      C = myAtom("C",i-1,aa_name_convert(string(prev)));
      N = myAtom("N",i,aa_name_convert(string(rescode)));
      add!(g,C,N);
   end

   # protein graph is ready
   return g;
end

# enriching the protein graph with edges connecting atoms separated by max 2 bonds
function enrich!(g::UndirectedGraph{myAtom})
   newedges = Vector{Tuple{myAtom,myAtom}}();

   # looking for the new edges to include
   for a1 in g.V
      for a2 in neighbors(g,a1)
         # edges (a1,a2) are, by definition, already present
         for a3 in neighbors(g,a2)
            if a3 != a1
               # adding the edge (a1,a3), if not present yet
               if !((a1,a3) in g.E) && !((a3,a1) in g.E)
                  push!(newedges,(a1,a3));
               end
               for a4 in neighbors(g,a3)
                  if a4 != a2 && a4 != a1
                     # adding the edge (a1,a4), if not present yet
                     if !((a1,a4) in g.E) && !((a4,a1) in g.E)
                        push!(newedges,(a1,a4));
                     end
                  end
               end
            end
         end
      end
   end

   # adding the new edges in the graph
   for e in newedges
      add!(g,e[1],e[2]);
   end

   # we are done
   return nothing;
end

