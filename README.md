# Monolis

- Monolithic domain decomposition based linear solver
- version 0.0.0

## Manual

https://morita.gitlab.io/monolis/

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
- Incomplete LU
- Jacobi
- Successive Over Relaxation

## Matrix storage format

- CSR format
    - n\*n blocking

## Interface

- Original

## Miscellanies

- MPI parallelization
- OpenMP parallelization
- Reordering with METIS version 5

### Compile flags

```
./install_lib.sh
make FLAGS=MPI,OPENMP,METIS,MUMPS
```

- MPI: Enable MPI parallelization
- OPENMP: Enable OpenMP parallelization
- METIS: Enable METIS
- MUMPS: Enable MUMPS
- DEBUG: Enable DEBUG

## License

- MIT

## Acknowledgements

- I am grateful to Mr. Yu IHARA for great assistance with the implementation of the fill-in determinaion in direct method.
- I would like to thank FrontISTR commons and FrontISTR for useful supports.

