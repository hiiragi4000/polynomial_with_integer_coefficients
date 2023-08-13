# Polynomial with Integer Coefficients

**Warning: The functions in this repository are UNSAFE IN ANY CONCURRENT SITUATION and thus are OBSOLETE now.
One should use _polynomial.hh_ in [this repository](https://github.com/hiiragi4000/cp_algos) instead.**

## Introduction

This library gives an implementation of the additions, subtractions, multiplications, and divisions of any two polynomials with integer-type coefficients.
One can pass an extra argument `m` to a constructor to indicate that the coefficients are in ring $\mathbb{Z}\_m$ instead of $\mathbb{Z}$.

To lower the time complexity of multiplications (as well as divisions) into $O(n\log n)$, we introduce **fast Fourier transform** (FFT) and **number-theoretic transform** (NTT).
The prime $p = 15\times2^{27}+1$ is chosen in NTT and is defined as the constant `nttp` in _polynomial.h_.
