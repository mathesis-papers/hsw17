# The weighted Ev character-average identity

Haobo Ma (ChronoAI Pte Ltd, Singapore) and Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore).

The note proves Hemmer–Straub–Westrem's Conjecture 17, equation (6.1), page 15 of *Electronic Journal of Combinatorics* 33(1) (2026), P1.57, DOI 10.37236/14401. Both centralizer-weighted sums vanish for odd positive n and equal binomial(n/2+2, 2) for even n. The pointwise identity also implies the N=1 case of Amdeberhan's q-series Conjecture 18.

## Files and build

- `note.tex`: complete LaTeX source, with bibliography included.
- `metadata.json`: scientific description, creators, license and source links.
- `verify.py`: dependency-free exact finite checks for Python 3.10 or later.
- `verification.txt`: the checker's concise expected output.

Build with Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, run a standard LaTeX installation twice; no BibTeX step is needed:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

Run the finite checks:

```sh
python3 verify.py > verification.actual.txt
diff -u verification.txt verification.actual.txt
```

## Scope

The written proof uses the Frobenius character formula, Patel's three-variable Schur coefficient detector and pointwise identity, and a formal cycle-index generating function. The all-degree weighted result is an explicit Ev-character application of the two-binary-quadratic invariant Hilbert series evaluated by Benanti et al. (2012), Example 3.3, journal page 160. The note gives its mathematical identification and an alternative coefficient calculation, not a first evaluation of that Hilbert series.

The declaration `result` in each of the following pinned files proves the corresponding coefficient-model statement:

- [HemmerStraubWestremWeightedAverage.lean](https://github.com/the-omega-institute/trureturing/blob/439c57ca26fd0c0d2de1f48bc492717bc6f7597c/D5/S3/Combinatorics/Partitions/HemmerStraubWestremWeightedAverage.lean): both weighted equalities for n >= 1.
- [HemmerStraubWestrem.lean](https://github.com/the-omega-institute/trureturing/blob/439c57ca26fd0c0d2de1f48bc492717bc6f7597c/D5/S3/Combinatorics/Partitions/HemmerStraubWestrem.lean): the pointwise identity for every positive-size partition.

In these sources `frobChar` is the coefficient of the Vandermonde times the power-sum product. Its interpretation as an ordinary irreducible character is the classical Frobenius formula; an internal Specht-module interpretation is separate from these declarations. The q-series corollary is a written deduction.

The finite checks compare Murnaghan–Nakayama character sums, the Laurent expression and the weighted values on all 67 partitions with 0 <= n <= 8, comprising 2101 distinct shape/cycle-type pairs. Selected identity-class values also match the hook-length formula. A second check compares every Laurent coefficient of the cycle-index recurrence with the three geometric factors through degree 60. The q-series coefficients are checked through degree 10. Arithmetic is integer or exact rational throughout. These finite ranges are distinct from the all-degree proofs.

## Consequences and attribution

The pointwise identity and same-parity Schur detector are due to Shivam Patel's public MathDB solution of problem 370115, Solution 1, Sections 1–3 (2026). The note states and proves the coefficient argument in full and derives the weighted consequences.

The Riordan/tableau formula follows already from the invariant-tensor decomposition of Meilhan–Suzuki (2015), Lemma 3.7, arXiv:1507.04454v1, and also from Okada (2016), Corollary 5.5(3), DOI 10.37236/6214. Jagenteufel (2018), DOI 10.37236/7713, Section 3.2.1 and Lemma 36, later gave a bijective proof; Hemmer–Straub–Westrem (2025) gave an independent enumerative proof. Callan's adjacent-trinomial formula is recorded in his OEIS comment dated 2004 and proved in his 2006 manuscript; Wallach's OEIS comment dated 2014 gives the rectangular all-ones case. The note rederives the single-part-size formulas, Westrem's first-row vanishing, and complete written consequences for mixed common divisors and arbitrary multiplicative weights. It makes no priority claim for these deductions.

The full Motzkin-level bijection question is answered affirmatively by Amaral, Gil and Weiner (22 September 2026), arXiv:2609.27096v1, Theorem 4.1. Their full Motzkin/tableau bijection restricts to Riordan paths and equal-parity tableaux. This cited result is distinct from Jagenteufel's earlier Riordan-only correspondence.

Questions 4 and 5 have affirmative answers for N=1 and N=2. The N=2 written proof identifies the adjoint representations of SO(5) and Sp(4), using an explicit symplectic exterior-square construction and the orthogonal/symplectic invariant criteria of Drensky–Hristova (2017), Corollaries 3.3 and 2.3. It also proves the N=2 q-series identity. These are written proofs, distinct from the two linked N=1 Lean declarations. The unrestricted equalities already have negative answers in the posing source at N=3; the checker reproduces its pointwise Appendix A example. The note makes no conclusion about the remaining N>=3 classifications or the q-series identity for N>=3.

Additional exact finite checks cover every monomial-symmetric basis element through degree 18 for the same-parity detector; Riordan-path counts against hook-length tableau degrees through degree 30; equal-part formulas for repetitions 1 through 15 and part values 1 through 5; every applicable common-divisor partition of sizes 1 through 12; every first-row vanishing pair through n=5; and the signed multiplicative weights w_k=(-1)^k(k+1) through degree 10. All ranges are finite illustrations of the separate written proofs.
The N=2 pointwise identity is checked on all 30 partitions with 0<=n<=6 by direct Murnaghan–Nakayama character sums.

## License

Copyright 2026 Haobo Ma and Wenlin Zhang. Licensed under Creative Commons Attribution 4.0 International (CC BY 4.0), identifier `cc-by-4.0`.
