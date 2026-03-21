import sympy as sp

# 1. System parameters for gap 3 base case (N=9)
n = 9
gap = 3
T = [i for i in range(n) if i % gap != 0] # [1, 2, 4, 5, 7, 8]
dim = len(T) # 6

def discrete_delta(i, j):
    return sp.Rational(1 if i % n == j % n else 0)

def B_coeff(k, j, l):
    # L96 Transverse Jacobian modular arithmetic
    d1 = discrete_delta(j, k-1) * discrete_delta(l, k-2)
    d2 = discrete_delta(j, k+2) * discrete_delta(l, k+1)
    d3 = discrete_delta(j, k+1) * discrete_delta(l, k+2)
    d4 = discrete_delta(j, k+1) * discrete_delta(l, k-1)
    return d1 - d2 + d3 - d4

# Generate the original M_k matrices (6x6)
M_matrices = []
for k in [3, 6, 9]:
    M = sp.zeros(dim, dim)
    for r, j in enumerate(T):
        for c, l in enumerate(T):
            M[r, c] = B_coeff(k, j, l)
    M_matrices.append(M)

# 2. Define the null vector w_0
# For T = [1, 2, 4, 5, 7, 8], the 3j+2 indices are 1, 3, 5.
w_0 = sp.Matrix([0, 1, 0, 1, 0, 1])

def compute_lie_algebra_rank(base_matrices, vector_length, max_depth):
    def compute_bracket(A, B):
        return A * B - B * A

    current_matrices = base_matrices.copy()
    all_matrices = base_matrices.copy()

    for depth in range(2, max_depth + 1):
        print(f"  Computing depth {depth}...") 
        new_matrices = []
        for A in base_matrices:
            for B in current_matrices:
                C = compute_bracket(A, B)
                if not C.is_zero_matrix:
                    new_matrices.append(C)
        all_matrices.extend(new_matrices)
        current_matrices = new_matrices
        print(f"    -> Generated {len(new_matrices)} new matrices ({len(all_matrices)} total)")
    
    print("  Stacking matrices to compute rank (this may take a few seconds)...")
    vectors = [M.reshape(vector_length, 1) for M in all_matrices]
    Lie_algebra_matrix = sp.Matrix.hstack(*vectors)
    return Lie_algebra_matrix.rank()

print("="*60)
print("PART 1: Full Space Lie Algebra (6x6)")
print("="*60)
print(f"Transverse Dimension: {dim}")
print(f"Target Rank for sl(6): {dim**2 - 1} (35)")

full_rank = compute_lie_algebra_rank(M_matrices, 36, max_depth=6)
print(f"Computed Rank on full H_I^perp: {full_rank}")
if full_rank < 35:
    print(f"CONCLUSION: The Lie algebra falls short of full rank by {35 - full_rank} dimensions.")
    print("This is due to the presence of the invariant subspace span(w_0).\n")


print("="*60)
print("PART 2: Quotient Space Lie Algebra (5x5)")
print("="*60)
# 3. Create a change of basis to isolate the quotient space
# We construct a basis of 5 linearly independent vectors, plus w_0
basis = [
    sp.Matrix([1, 0, 0, 0, 0, 0]),  # e_1
    sp.Matrix([0, 0, 1, 0, 0, 0]),  # e_4
    sp.Matrix([0, 0, 0, 0, 1, 0]),  # e_7
    sp.Matrix([0, 1, 0, -1, 0, 0]), # 3j+2 difference
    sp.Matrix([0, 0, 0, 1, 0, -1]), # 3j+2 difference
    w_0                             # The null vector
]
P = sp.Matrix.hstack(*basis)
P_inv = P.inv()

# 4. Project M_k onto the 5D quotient space
M_tilde = []
for M in M_matrices:
    M_transformed = P_inv * M * P
    # The top-left 5x5 block is the induced map on the quotient space
    M_tilde.append(M_transformed[:5, :5])

print(f"Quotient Space Dimension: 5")
print(f"Target Rank for sl(5): {5**2 - 1} (24)")

quotient_rank = compute_lie_algebra_rank(M_tilde, 25, max_depth=5)
print(f"Final Computed Rank on quotient space: {quotient_rank}")
if quotient_rank == 24:
    print("SUCCESS: The matrices span exactly sl(5) on the quotient space!")
else:
    print("FAILURE: The matrices do not span the full quotient space algebra.")
