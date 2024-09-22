
# Distance Geometry repository

This repository is meant to contain several software projects aiming at studying
and solving the Distance Geometry Problem (DGP).

## Subset Sum Problem (SSP)

This is a C implementation of the Branch-and-Prune (BP) algorithm, initially developed
for the solution of [discretizable DGPs](https://link.springer.com/article/10.1007/s11590-011-0358-3),
which is here adapted to solve instances of the **Subset Sum Problem** (SSP). The SSP
is one of the very well-known classical problems in combinatorial optimization, which
belong to the **NP**-complete complexity class. Even if it has an exponential worst-case
complexity, BP is able to find solutions to the hardest instances of the SSP (those for
which the *density* index is very close to 1, see function ```ssp_gen_hard```) in a 
reasonable amount of time. The implementation of other algorithms, claiming a *pseudo*
polynomial complexity on simpler SSP instances, is currently under study. The implementation
of BP which is already on the repository (see this [C file](./ssp/bp.c) is actually quite
difficult to beat, because it uses no extra memory, and the local variables necessary for 
the recursive calls are allocated directly on the "stack", hence implying a very fast
memory access.

## juliaDGP

This is an initial implementation of a generic software tool for the DGP, written
in [Julia](https://julialang.org/) programming language. The current preliminary 
version is basically the result of conception work, where the well-known object-oriented 
paradigm is rather replaced by the *multiple-dispatch* paradigm exploited in Julia.
The initial focus is on the typical DGP application in structural biology: this
tool is "aware" of the structure of all amino acids that can be involved in the
protein synthesis. This current version is already capable of creating DGP instances 
from PDB files, as well as from STAR files, and to perform some basic comparisons 
through the use of the newly introduced ```ErrorList``` type. Some documentation 
for the use of the software tool will be added to the repository in the near future.

Meanwhile, in order to load the source code in Julia and explore it, just run 
(in Julia):

	julia> include("JuliaDGP.jl")

and read the one-line comments at the top of every function before invoking them.

## javaCMP

This a Java project (named ```javaCMP```) that implements a *multi-representation* 
system for a given simple directed graph, which is at the basis of the Coherent
Multirepresentation Problem (CMP). The CMP is a possible extension of the DGP, 
which was introduced in this [conference paper](https://link.springer.com/chapter/10.1007/978-3-031-34953-9_27).

------------------------------------------------------------------------------------

