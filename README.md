# Monolis

- Morita's non-overlapping / overlapping domain decomposition based linear equation solver

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
- Reordering with metis version 5

## Compile

Confirm and set the following variables in Makefile.

```
FC     = mpif90
FFLAGS = -O2
METIS_DIR = $(path to metis)
METIS_INC = -I $(METIS_DIR)/include
METIS_LIB = -L$(METIS_DIR)/lib -lmetis
```

### with MPI and METIS

```
make FLAGS=MPI,METIS
```

### with MPI

```
make FLAGS=MPI
```

### with METIS

```
make FLAGS=METIS
```

## License

- MIT

## Acknowledgements

- I am grateful to Mr. Yu IHARA for great assistance with the implementation of the fill-in determinaion in direct method.
- I would like to thank FrontISTR commons and FrontISTR for useful supports.

