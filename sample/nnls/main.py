import numpy as np
from scipy.linalg import solve, LinAlgWarning
import warnings

def nnls(A, b, maxiter=None, *, atol=None):
    A = np.asarray_chkfinite(A)
    b = np.asarray_chkfinite(b)

    if len(A.shape) != 2:
        raise ValueError("Expected a two-dimensional array (matrix)" +
                         f", but the shape of A is {A.shape}")
    if len(b.shape) != 1:
        raise ValueError("Expected a one-dimensional array (vector)" +
                         f", but the shape of b is {b.shape}")

    m, n = A.shape

    if m != b.shape[0]:
        raise ValueError(
                "Incompatible dimensions. The first dimension of " +
                f"A is {m}, while the shape of b is {(b.shape[0], )}")

    x, rnorm, mode = _nnls(A, b, maxiter, tol=atol)
    if mode != 1:
        raise RuntimeError("Maximum number of iterations reached.")

    return x, rnorm


def _nnls(A, b, maxiter=None, tol=None):
    """
    This is a single RHS algorithm from ref [2] above. For multiple RHS
    support, the algorithm is given in  :doi:`10.1002/cem.889`
    """
    m, n = A.shape

    AtA = A.T @ A
    Atb = b @ A  # Result is 1D - let NumPy figure it out

    if not maxiter:
        maxiter = 3*n
    if tol is None:
        tol = 10 * max(m, n) * np.spacing(1.)

    # Initialize vars
    x = np.zeros(n, dtype=np.float64)
    s = np.zeros(n, dtype=np.float64)
    # Inactive constraint switches
    P = np.zeros(n, dtype=bool)

    # Projected residual
    w = Atb.copy().astype(np.float64)  # x=0. Skip (-AtA @ x) term

    # Overall iteration counter
    # Outer loop is not counted, inner iter is counted across outer spins
    iter = 0

    while (not P.all()) and (w[~P] > tol).any():  # B
        # Get the "most" active coeff index and move to inactive set
        k = np.argmax(w * (~P))  # B.2
        P[k] = True  # B.3

        # Iteration solution
        s[:] = 0.
        # B.4
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='Ill-conditioned matrix',
                                    category=LinAlgWarning)
            s[P] = solve(AtA[np.ix_(P, P)], Atb[P], assume_a='sym', check_finite=False)

        # Inner loop
        while (iter < maxiter) and (s[P].min() < 0):  # C.1
            iter += 1
            inds = P * (s < 0)
            alpha = (x[inds] / (x[inds] - s[inds])).min()  # C.2
            x *= (1 - alpha)
            x += alpha*s
            P[x <= tol] = False
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', message='Ill-conditioned matrix',
                                        category=LinAlgWarning)
                s[P] = solve(AtA[np.ix_(P, P)], Atb[P], assume_a='sym',
                             check_finite=False)
            s[~P] = 0  # C.6

        x[:] = s[:]
        w[:] = Atb - AtA @ x

        if iter == maxiter:
            # Typically following line should return
            # return x, np.linalg.norm(A@x - b), -1
            # however at the top level, -1 raises an exception wasting norm
            # Instead return dummy number 0.
            return x, 0., -1

    return x, np.linalg.norm(A@x - b), 1

A = np.array([[1, 0, 1, 0, 0], 
              [1, 2, 0, 0, 0], 
              [2, 0, 0, 1, 0], 
              [1, 0, 4, 0, 0], 
              [3, 0, 0, 0, 1], 
              [1, 3, 0, 0, 0], 
              [4, 0, 1, 2, 0], 
              [0, 1, 0, 1, 1]])

b = np.array([-1,-1, 2, 2,-3, -3, 4, -4])

x = nnls(A, b)

print(x)

