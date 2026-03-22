# Supplementary Scripts for Lorenz 96 Lie Algebra Generation

This repository contains supplementary verification scripts for the paper:

**"Non-uniqueness of stationary measures for stochastic systems with almost surely invariant manifolds"**
by Jacob Bedrossian, Alex Blumenthal, and Sam Punshon-Smith.

## Overview

The L96 model is forced on every 4th mode ($I = 4\mathbb{Z}/N\mathbb{Z}$, gap 4). The Lie algebra generation result — that $\mathrm{Lie}(\{M_k : k \in I\}) = \mathfrak{sl}(H_I^\perp)$ for all $N = 4K$, $K \ge 3$ — is proved **algebraically** in Appendix C of the paper, exploiting a tensor decomposition $H_I^\perp \simeq \mathbb{R}^K \otimes \mathbb{R}^3$. The scripts in this repository provide:

1. **Independent numerical verification** of every bracket computation in the algebraic proof.
2. **Symbolic verification** of the gap-3 quotient space analysis (Remark 4.2).
3. **Symbolic eigenvalue computation** for the transverse linearization instability (Section 5).

## Files

| File | Description |
|---|---|
| `verify_algebraic_proof.py` | Numerically verifies all 77 bracket computations in the algebraic proof (Appendix C). Checks Lemma C.2 (the seed), Lemma C.3 (macroscopic shifts), Lemma C.4 (local algebra generation), and the final proposition assembly. Uses `numpy`. |
| `gap3_quotient_bracket.py` | Symbolically verifies the gap-3 case (Remark 4.2). Computes the Lie bracket rank on both the full 6D transverse space (rank 29/35, showing obstruction) and the 5D quotient space (rank 24/24, confirming generation). Uses `sympy`. |
| `symbolic_eigen.py` | Symbolically computes the $9 \times 9$ matrix $DB(y)|_{H_I^\perp}$ for $y = a(e_0 + e_4 + e_8)$ and factors its characteristic polynomial, verifying the instability used in Lemma 5.7. Uses `sympy`. |
| `old/` | Contains the original computer-assisted proof notebook from earlier versions of the paper. |

## Requirements

- Python 3
- `numpy`
- `sympy`

## Usage

```bash
# Verify every computation in the algebraic proof (Appendix C)
python3 verify_algebraic_proof.py

# Verify gap-3 full space and quotient space bracket span (Remark 4.2)
python3 gap3_quotient_bracket.py

# Verify characteristic polynomial factorization (Section 5)
python3 symbolic_eigen.py
```