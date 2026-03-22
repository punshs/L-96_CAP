"""
Comprehensive verification of every computation in the proof of app_cap.tex.
Uses 0-indexed internal matrices but reports in 1-indexed notation.
"""
import numpy as np
import sys

def E(n,i,j):
    M=np.zeros((n,n)); M[i,j]=1.0; return M
def kron(A,B): return np.kron(A,B)
def bracket(A,B): return A@B - B@A
def mat_eq(A,B,tol=1e-10): return np.max(np.abs(A-B))<tol

K=3; n=3*K
# Internal matrices (0-indexed)
A3=E(3,2,1); B3=E(3,0,1)-E(3,1,0); C3=-E(3,0,2)
E31=E(3,2,0); E11=E(3,0,0); E32=E(3,2,1); E12=E(3,0,1)
E21=E(3,1,0); E13=E(3,0,2); E33=E(3,2,2)
E22=E(3,1,1)

def Mk(m):
    return kron(E(K,(m-1)%K,(m-1)%K),A3)+kron(E(K,m%K,m%K),B3)+kron(E(K,m%K,(m-1)%K),C3)

ok_count=0; fail_count=0
def check(name, result, expected):
    global ok_count, fail_count
    if mat_eq(result, expected):
        print(f"  ✓ {name}")
        ok_count += 1
    else:
        print(f"  ✗ FAIL: {name}")
        fail_count += 1
        # Show difference
        diff = result - expected
        print(f"    Max diff: {np.max(np.abs(diff))}")

print("="*70)
print("LEMMA C.2: The Seed")
print("="*70)

# Line 71: AC = 0
check("AC = 0", A3@C3, np.zeros((3,3)))

# Line 71: CB = 0
check("CB = 0", C3@B3, np.zeros((3,3)))

# Line 71: C^2 = 0
check("C^2 = 0", C3@C3, np.zeros((3,3)))

# Line 77: [B,A] = E_{3,1}
check("[B,A] = E_{3,1}", bracket(B3, A3), E31)

# Line 77: Check the intermediate: (E12-E21)*E32 = 0
check("(E12-E21)*E32 = 0", B3@A3, np.zeros((3,3)))

# Line 77: E32*(E12-E21) = -E31
check("E32*(E12-E21) = -E31", A3@B3, -E31)

# Line 62: Full bracket [M_4m, M_4(m+1)] = E_{m,m} x E_{3,1}
for m in range(K):
    expected = kron(E(K,m,m), E31)
    check(f"[M_{m}, M_{(m+1)%K}] = E_{{{m},{m}}}⊗E_31", bracket(Mk(m), Mk((m+1)%K)), expected)

print()
print("="*70)
print("LEMMA C.3: Macroscopic Shifts")
print("="*70)

# Line 93: [E31, A] = [E31, E32] = 0
check("[E31, E32] = 0", bracket(E31, A3), np.zeros((3,3)))

# Line 93: E31*E32 = 0
check("E31*E32 = 0", E31@A3, np.zeros((3,3)))

# Line 93: E32*E31 = 0
check("E32*E31 = 0", A3@E31, np.zeros((3,3)))

# Line 98: C*E31 = (-E13)(E31) = -E11
check("C*E31 = -E11", C3@E31, -E11)

# Line 100: Full bracket [X_{m-1}, M_{4m}] = Z_m = E_{m,m-1} x E11
for m in range(K):
    Xm1 = kron(E(K,(m-1)%K,(m-1)%K), E31)
    expected = kron(E(K,m,(m-1)%K), E11)
    check(f"[X_{(m-1)%K}, M_{m}] = E_{{{m},{(m-1)%K}}}⊗E_11", bracket(Xm1, Mk(m)), expected)

# Line 105: [Z_l, Z_{l-1}] = E_{l,l-2} x E11
for l in range(K):
    Zl = kron(E(K,l,(l-1)%K), E11)
    Zl1 = kron(E(K,(l-1)%K,(l-2)%K), E11)
    expected = kron(E(K,l,(l-2)%K), E11)
    check(f"[Z_{l}, Z_{(l-1)%K}] = E_{{{l},{(l-2)%K}}}⊗E_11", bracket(Zl, Zl1), expected)

# Line 113: [E_{m,m+1}xE11, E_{m+1,m}xE11] = (E_{mm}-E_{m+1,m+1})xE11
for m in range(K):
    fwd = kron(E(K,m,(m+1)%K), E11)
    bwd = kron(E(K,(m+1)%K,m), E11)
    expected = kron(E(K,m,m)-E(K,(m+1)%K,(m+1)%K), E11)
    check(f"[E_{{{m},{(m+1)%K}}}⊗E11, E_{{{(m+1)%K},{m}}}⊗E11] = Cartan", bracket(fwd, bwd), expected)

print()
print("="*70)
print("LEMMA C.4: Local Algebra Shattering")
print("="*70)

# Line 124: [E31, B] = E32
check("[E31, B] = [E31, E12-E21] = E32", bracket(E31, B3), E32)

# Line 124: E31*C = E31*(-E13) = -E33
check("E31*C = -E33", E31@C3, -E33)

# Line 124: Full [X_m, M_4m]
for m in range(K):
    Xm = kron(E(K,m,m), E31)
    expected = kron(E(K,m,m), E32) - kron(E(K,m,(m-1)%K), E33)
    check(f"[X_{m}, M_{m}] = E_{{{m},{m}}}⊗E32 - E_{{{m},{(m-1)%K}}}⊗E33", bracket(Xm, Mk(m)), expected)

# Line 128: ad^2_{M4m}(Xm) = -E_{m,m}xE31 - E_{m,m-1}xE32
for m in range(K):
    Xm = kron(E(K,m,m), E31)
    inner = bracket(Mk(m), Xm)
    ad2 = bracket(Mk(m), inner)
    expected = -kron(E(K,m,m), E31) - kron(E(K,m,(m-1)%K), E32)
    check(f"ad^2_M{m}(X_{m}) = -E_{{{m},{m}}}⊗E31 - E_{{{m},{(m-1)%K}}}⊗E32", ad2, expected)

# Line 132: X_m + ad^2 = -E_{m,m-1}xE32
for m in range(K):
    Xm = kron(E(K,m,m), E31)
    inner = bracket(Mk(m), Xm)
    ad2 = bracket(Mk(m), inner)
    result = Xm + ad2
    expected = -kron(E(K,m,(m-1)%K), E32)
    check(f"X_{m} + ad^2 = -E_{{{m},{(m-1)%K}}}⊗E32", result, expected)

# Line 138: CE32 = (-E13)E32 = -E12
check("CE32 = (-E13)E32 = -E12", C3@E32, -E12)

# Line 138: E32^2 = 0
check("E32^2 = 0", E32@E32, np.zeros((3,3)))

# Line 136: [E_{m,m-1}xE32, M_{4(m+1)}] = -E_{m+1,m-1}xCE32 = E_{m+1,m-1}xE12
for m in range(K):
    elem = kron(E(K,m,(m-1)%K), E32)
    expected = kron(E(K,(m+1)%K,(m-1)%K), E12)
    check(f"[E_{{{m},{(m-1)%K}}}⊗E32, M_{(m+1)%K}] = E_{{{(m+1)%K},{(m-1)%K}}}⊗E12", bracket(elem, Mk((m+1)%K)), expected)

# Line 144: [E_{m-1,m+1}xE11, E_{m+1,m-1}xE12] = E_{m-1,m-1}xE12
for m in range(K):
    shift = kron(E(K,(m-1)%K,(m+1)%K), E11)
    elem = kron(E(K,(m+1)%K,(m-1)%K), E12)
    expected = kron(E(K,(m-1)%K,(m-1)%K), E12)
    check(f"shift E_{{{(m+1)%K},{(m-1)%K}}}⊗E12 to diagonal", bracket(shift, elem), expected)

# Line 144: E11*E12 = E12
check("E11*E12 = E12", E11@E12, E12)

# Line 148: [E_{m,m}xE31, E_{m,m}xE12] = E_{m,m}xE32
check("[E31, E12] = E32", E31@E12, E32)

# Line 150: Mtilde = M - E_{m-1,m-1}xE32 = E_{m,m}xB + E_{m,m-1}xC
for m in range(K):
    Mtilde = Mk(m) - kron(E(K,(m-1)%K,(m-1)%K), E32)
    expected = kron(E(K,m,m), B3) + kron(E(K,m,(m-1)%K), C3)
    check(f"Mtilde_{m} = E_{{{m},{m}}}⊗B + E_{{{m},{(m-1)%K}}}⊗C", Mtilde, expected)

# Line 152: E12*(-E13) = 0 (C-block annihilation)
check("E12*(-E13) = 0", E12@C3, np.zeros((3,3)))

# Line 154: [B, E12] = [E12-E21, E12] = E11-E22
check("[E12-E21, E12] = E11-E22", bracket(B3, E12), E11-E22)

# Line 154: Full H_m = [Mtilde, E_{m,m}xE12]
for m in range(K):
    Mtilde = Mk(m) - kron(E(K,(m-1)%K,(m-1)%K), E32)
    elem = kron(E(K,m,m), E12)
    expected = kron(E(K,m,m), E11-E22)
    check(f"H_{m} = [Mtilde_{m}, E_{{{m},{m}}}⊗E12] = E_{{{m},{m}}}⊗(E11-E22)", bracket(Mtilde, elem), expected)

# Line 156: [E11-E22, E12-E21] = 2E12+2E21
H_int = E11-E22
check("[E11-E22, E12-E21] = 2E12+2E21", bracket(H_int, B3), 2*E12+2*E21)

# Line 156: (E11-E22)(-E13) = -E13
check("(E11-E22)(-E13) = -E13", H_int@C3, C3)

# Line 158: [H_m, Mtilde] full check
for m in range(K):
    Mtilde = Mk(m) - kron(E(K,(m-1)%K,(m-1)%K), E32)
    Hm = kron(E(K,m,m), H_int)
    expected = kron(E(K,m,m), 2*E12+2*E21) + kron(E(K,m,(m-1)%K), C3)
    check(f"[H_{m}, Mtilde_{m}] = E_{{{m},{m}}}⊗(2E12+2E21) + E_{{{m},{(m-1)%K}}}⊗(-E13)", bracket(Hm, Mtilde), expected)

# Line 160: [H_m, Mtilde] - Mtilde = E_{m,m}x(E12+3E21)
for m in range(K):
    Mtilde = Mk(m) - kron(E(K,(m-1)%K,(m-1)%K), E32)
    Hm = kron(E(K,m,m), H_int)
    result = bracket(Hm, Mtilde) - Mtilde
    expected = kron(E(K,m,m), E12+3*E21)
    check(f"[H_{m}, Mtilde_{m}] - Mtilde_{m} = E_{{{m},{m}}}⊗(E12+3E21)", result, expected)

# Line 164: Mtilde - E_{m,m}x(E12-E21) = -E_{m,m-1}xE13
for m in range(K):
    Mtilde = Mk(m) - kron(E(K,(m-1)%K,(m-1)%K), E32)
    result = Mtilde - kron(E(K,m,m), B3)
    expected = kron(E(K,m,(m-1)%K), C3)  # = -E_{m,m-1}xE13
    check(f"Mtilde_{m} - E_{{{m},{m}}}⊗B = -E_{{{m},{(m-1)%K}}}⊗E13", result, expected)

# Line 166: Shift -E_{m,m-1}xE13 to diagonal
for m in range(K):
    bwd = kron(E(K,(m-1)%K,m), E11)
    elem = -kron(E(K,m,(m-1)%K), E13)
    result = bracket(bwd, elem)
    expected = -kron(E(K,(m-1)%K,(m-1)%K), E13)
    check(f"[E_{{{(m-1)%K},{m}}}⊗E11, -E_{{{m},{(m-1)%K}}}⊗E13] = -E_{{{(m-1)%K},{(m-1)%K}}}⊗E13", result, expected)

print()
print("="*70)
print("PROOF OF PROPOSITION (Final Assembly)")
print("="*70)

# Line 174: [E_{a,b}xE11, E_{b,b}xE_{1,beta}] = E_{a,b}xE_{1,beta}
for beta in range(3):
    a,b = 0,1
    shift = kron(E(K,a,b), E11)
    local = kron(E(K,b,b), E(3,0,beta))
    expected = kron(E(K,a,b), E(3,0,beta))
    check(f"[E_{{{a},{b}}}⊗E11, E_{{{b},{b}}}⊗E_{{1,{beta+1}}}] = E_{{{a},{b}}}⊗E_{{1,{beta+1}}}", bracket(shift, local), expected)

# Line 178: [E_{a,a}xE_{alpha,1}, E_{a,b}xE_{1,beta}] = E_{a,b}xE_{alpha,beta}
for alpha in range(3):
    for beta in range(3):
        a,b = 0,1
        local = kron(E(K,a,a), E(3,alpha,0))
        elem = kron(E(K,a,b), E(3,0,beta))
        expected = kron(E(K,a,b), E(3,alpha,beta))
        check(f"[E_{{{a},{a}}}⊗E_{{{alpha+1},1}}, E_{{{a},{b}}}⊗E_{{1,{beta+1}}}] = E_{{{a},{b}}}⊗E_{{{alpha+1},{beta+1}}}", bracket(local, elem), expected)

print()
print("="*70)
print(f"SUMMARY: {ok_count} passed, {fail_count} failed")
print("="*70)
sys.stdout.flush()
