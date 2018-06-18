# Monolis

- Morita's non-overlapping / overlapping domain decomposition based linear equation solver

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

## License

- MIT
