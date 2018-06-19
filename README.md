# Monolis

- Morita's non-overlapping / overlapping domain decomposition based linear equation solver

## Direct solver

- LU factorization (LU)

## Iterative solver

- Conjugate Gradient Method (CG)
- Biconjugate Gradient Stabilized Method (BiCGSTAB)
- Biconjugate Gradient Stabilized Method without preconditioning (BiCGSTAB-noprec)
- Grop's Conjugate Gradient Method (GropCG)
- Pipelined Conjugate Gradient Method (PipeCG)
- Pipelined Conjugate Residual Method (PipeCR)
- Pipelined Biconjugate Gradient Stabilized Method (PipeBiCGSTAB)
- Pipelined Biconjugate Gradient Stabilized Method without preconditioning (PipeBiCGSTAB-noprec)
- Communication-avoiding Biconjugate Gradient Stabilized Method without preconditioning (CABiCGSTAB-noprec)
- Successive Over Relaxation (SOR)
- Iterative Refinement (IR)

## Preconditioning

- Diagonal scaling

## Matrix storage format

- CSR format
    - n\*n blocking
    - storing a sparse matrix separately as A = L + D + U

## Interface

- Original

## Miscellanies

- MPI parallelization (-DWITH_MPI)
- Reordering with metis version 5 (-DWITH_METIS)

## License

- MIT

## Acknowledgements

- I am grateful to Mr. Yu IHARA for great assistance with the implementation of the fill-in determinaion in direct method.
- I would like to thank FrontISTR commons and FrontISTR for useful supports.

