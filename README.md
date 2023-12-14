# This codebase consists of codes for the following:
## [1] Fortran and Python codes for generating hamiltonian matrix for various few/many body models (including long range hopping, interactions, et.c.)
## [2] Graph Traversal Algorithms using mathematica and Python package - NetworkX.

## [3] Block Tri-Diagonalisation of Sparse, Banded Matrix.
The red lines represent the envelops of the sparse, banded matrices. Algorithm to derive a block tridiagonal representation from the envelop profile have been implemented.

<img width="450" alt="bandsize" src="https://github.com/karakashdeep/many-few-body-model-simulations/assets/60542977/fff447a2-5776-4e62-9ede-5d22bf52b0f4">

## [4] Implementation of [Fock space Recursive Greens functions (FRGF)](https://iopscience.iop.org/article/10.1088/1751-8121/acef7a/meta) and its integration with Bandwidth Reduction Algorithms
Recursive procedure involves obtaining a Block-Tridiagonal Representation, followed by (i) __Forward__ (ii) __Backward__ and (iii) __off-Diagonal__ Recursive relations to obtain the complete resolvent $[(\omega+i\eta)\mathcal{I}-\mathcal{H}]^{-1}$

<img width="600" alt="recursion" src="https://github.com/karakashdeep/many-few-body-model-simulations/assets/60542977/4e05b4bd-e9fd-4abe-9771-5a9096d2ea75">

## [5] Some prelim results and testing ideas.
### Memory Scaling (for 1d Spinless Hubbard Model)
The y-axis represents the _Dynamic Memory Requirement_ or the size of the largest matrix required to be stored at a time.

<img width="400" alt="image" src="https://github.com/karakashdeep/many-few-body-model-simulations/assets/60542977/0ca8314e-6418-4f5a-84e9-e38854fdd817">

### Minimum Spanning Tree with coresponding "Level sets" on the _4x4_ 2d lattice (TBM)
A way of Visualizing the scheme on a one-body problem i.e 2d Tight Binding Model. Figure on the left is a tree _spanning_ the _fock space graph_ of TBM. Figure on the right are "level set" contours on the 2d lattice

<img width="400" alt="Screenshot 2023-11-21 at 4 41 46 PM" src="https://github.com/karakashdeep/many-few-body-model-simulations/assets/60542977/440abb54-4004-4bf8-9c84-e217d36f4df9"> 

<img width="400" alt="Screenshot 2023-11-21 at 4 42 00 PM" src="https://github.com/karakashdeep/many-few-body-model-simulations/assets/60542977/5dad8648-68f0-42d9-b27e-f97c472f8193">

### __"Level sets"__ on the _6x6_ 2d lattice (TBM)
<img width="450" alt="Screenshot 2023-11-21 at 4 42 15 PM" src="https://github.com/karakashdeep/many-few-body-model-simulations/assets/60542977/fa3acc2b-c50c-4193-9856-cf526b323f71">

