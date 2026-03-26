# Supplementary Scripts for Lorenz 96 Lie Algebra Generation

This repository contains supplementary verification scripts for the paper:

**"Non-uniqueness of stationary measures for stochastic systems with almost surely invariant manifolds"**
by Jacob Bedrossian, Alex Blumenthal, and Sam Punshon-Smith.

## Overview

The L96 model is forced on every 4th mode ($I = 4\mathbb{Z}/N\mathbb{Z}$, gap 4). The Lie algebra generation result — that $\mathrm{Lie}(\{M_k : k \in I\}) = \mathfrak{sl}(H_I^\perp)$ for all $N = 4K$, $K \ge 3$ — is proved **algebraically** in Appendix C of the paper, exploiting a tensor decomposition $H_I^\perp \simeq \mathbb{R}^K \otimes \mathbb{R}^3$. No computer-assisted proof (CAP) is required for the Lie algebra generation argument.

The only computation in the paper that relies on computer assistance is the **factorization of the characteristic polynomial** in Lemma 5.7 (Section 5). The script `symbolic_eigen.py` performs this verification.

## Files

| File | Description |
|---|---|
| `symbolic_eigen.py` | Symbolically computes the $9 \times 9$ matrix $DB(y)\|_{H_I^\perp}$ for $y = a(e_0 + e_4 + e_8)$ and factors its characteristic polynomial, verifying the factorization stated in the proof of Lemma 5.7 (transverse linear instability). Uses `sympy`. |

## Requirements

- Python 3
- `sympy`

## Usage

```bash
# Verify characteristic polynomial factorization (Lemma 5.7)
python3 symbolic_eigen.py
```