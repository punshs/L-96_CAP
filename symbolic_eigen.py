"""
Symbolic verification of the characteristic polynomial of DB(y)|_{H_I^perp}
for the Lorenz-96 system with gap-4 forcing (N = 12, I = {0, 4, 8}).

This script accompanies Lemma 5.7 of:
  "Non-uniqueness of stationary measures for stochastic systems
   with almost surely invariant manifolds"
  by J. Bedrossian, A. Blumenthal, and S. Punshon-Smith.

For y = a(e_0 + e_4 + e_8), the transverse derivative DB(y)|_{H_I^perp}
reduces to a 9x9 block on the indices T_0 = {1,2,3,5,6,7,9,10,11}.
This script constructs that block symbolically and verifies the factored
characteristic polynomial stated in the proof.

Requirements: sympy
Usage: python symbolic_eigen.py
"""

from sympy import symbols, Matrix, det, factor, eye, expand

# ---------------------------------------------------------------------------
# 1. Build the 9x9 matrix M = DB(y)|_{H_I^perp} for y = e_0 + e_4 + e_8
# ---------------------------------------------------------------------------
#
# The Lorenz-96 nonlinearity is B(u, u)_j = u_{j+1} u_{j-1} - u_{j-2} u_{j-1},
# so the linearization acts on a perturbation w by
#
#     (DB(y) w)_j = (w_{j+1} - w_{j-2}) y_{j-1} + (y_{j+1} - y_{j-2}) w_{j-1}.
#
# With N = 12 and y_k = 1 if k in {0,4,8}, y_k = 0 otherwise, the matrix
# DB(y) restricted to H_I^perp is block-diagonal: a single 9x9 block on
# T_0 = {1,2,3,5,6,7,9,10,11} and zeros on the remaining transverse indices.
#
# We construct this block directly from the L96 formula.

N = 12
T = [j for j in range(N) if j % 4 != 0]     # transverse indices
idx = {j: i for i, j in enumerate(T)}         # map: mode index -> matrix row/col


def y(k):
    """Value of y_k for y = e_0 + e_4 + e_8 (indices mod N)."""
    return 1 if k % N in {0, 4, 8} else 0


M = Matrix.zeros(9, 9)

for row, j in enumerate(T):
    # Term 1: y_{j-1} * (w_{j+1} - w_{j-2})
    c1 = y(j - 1)
    if c1:
        jp1 = (j + 1) % N
        jm2 = (j - 2) % N
        if jp1 in idx:
            M[row, idx[jp1]] += c1
        if jm2 in idx:
            M[row, idx[jm2]] -= c1

    # Term 2: (y_{j+1} - y_{j-2}) * w_{j-1}
    c2 = y(j + 1) - y(j - 2)
    if c2:
        jm1 = (j - 1) % N
        if jm1 in idx:
            M[row, idx[jm1]] += c2

# ---------------------------------------------------------------------------
# 2. Display the matrix
# ---------------------------------------------------------------------------

print("Lorenz-96 gap-4 eigenvalue verification")
print("=" * 55)
print(f"\nN = {N},  I = {{0, 4, 8}},  T_0 = {T}")
print(f"\n9x9 matrix  M = DB(y)|_{{H_I^perp}}  for y = e_0 + e_4 + e_8:\n")

labels = [f"e_{j}" for j in T]
col_header = "       " + "  ".join(f"{s:>4}" for s in labels)
print(col_header)
for i, lbl in enumerate(labels):
    entries = "  ".join(f"{int(M[i, j]):4d}" for j in range(9))
    print(f"  {lbl:>4}: {entries}")

# ---------------------------------------------------------------------------
# 3. Compute and factor the characteristic polynomial of a*M
# ---------------------------------------------------------------------------

a, t = symbols('a t')

char_poly = det(a * M - t * eye(9))
char_poly_factored = factor(char_poly)

print(f"\nCharacteristic polynomial  det(aM - tI):")
print(f"  {char_poly_factored}")

# ---------------------------------------------------------------------------
# 4. Verify the claimed factorization from Lemma 5.7
# ---------------------------------------------------------------------------

claimed = -(a**2 - a*t + t**2) * (t**3 + a**2*t - a**3) * \
           (t**4 + a*t**3 + 2*a**2*t**2 + 2*a**3*t + a**4)

diff = expand(char_poly - claimed)

print(f"\nClaimed factorization from Lemma 5.7:")
print(f"  -(a^2 - at + t^2)(t^3 + a^2 t - a^3)(t^4 + at^3 + 2a^2 t^2 + 2a^3 t + a^4)")
print(f"\nVerification (difference = 0?):  {diff == 0}")

if diff == 0:
    print("\n*** Factorization VERIFIED. ***")
else:
    print(f"\n*** MISMATCH: difference = {diff} ***")
