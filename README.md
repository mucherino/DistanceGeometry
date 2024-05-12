
# Distance Geometry repository

This repository is meant to contain several software tools aiming at studying
and solving the Distance Geometry Problem (DGP).

## juliaDGP

This is an initial implementation of a generic software tool for the DGP, written
in [Julia](https://julialang.org/) programming language. The current preliminary 
version is basically the result of conception work, where the well-known object-oriented 
paradigm is rather replaced by the *multiple-dispatch* paradigm exploited in Julia.
The initial focus is on the typical DGP application in structural biology: this
tool is "aware" of the structure of all amino acids that can be involved in the
protein synthesis. This current version is already capable of creating DGP instances 
from PDB files, as well as from STAR files, and to perform some basic comparisons 
through the use of the newly introduced ```ErrorList``` type. A minimal documentation 
for the use of the software tool will be added to the repository in the near future.

Meanwhile, in order to load the source code in Julia and explore it, just run 
(in Julia):

	julia> include("JuliaDGP.jl")

and read the one-line comments at the top of every function before invoking them.

## javaCMP

This a Java project (named ```javaCMP```) that implements a *multi-representation* 
system for a given simple directed graph, which is at the basis of the Coherent
Multirepresentation Problem (CMP). The CMP is a possible extension of the DGP, 
which was introduced in the following article:

- A. Mucherino,
  *The Coherent Multi-Representation Problem with Applications in Structural Biology*,
  Lectures Notes in Computer Science 13919, Lecture Notes in Bioinformatics series, 
  Proceedings of IWBBIO23, Gran Canaria, Spain, 338-346, 2023.
  [SpringerLink](https://link.springer.com/chapter/10.1007/978-3-031-34953-9_27).
  The article is also available on HAL open archives.

------------------------------------------------------------------------------------

