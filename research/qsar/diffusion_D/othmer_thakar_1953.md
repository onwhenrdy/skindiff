---
name: othmer_thakar_1953
predicts: D
layer: aqueous
units: cm^2/s
descriptors: solute molar volume V_A, water viscosity eta_w
domain: dilute aqueous solution; small organic non-electrolytes
status: classical
---

# Othmer–Thakar 1953 — Empirical aqueous diffusion correlation

## Citation
Othmer DF, Thakar MS. "Correlating diffusion coefficients in liquids."
Ind Eng Chem. 1953;45(3):589–593. DOI: 10.1021/ie50519a036.

## Equation
Predates Wilke–Chang; the original water-only form:

```
D_aq [cm^2/s] = 14.0e-5 / ( eta_w^1.1 * V_A^0.6 )
```

with eta_w in cP, V_A in cm^3/mol (LeBas). At 25 C (eta_w = 0.89):

```
D_aq [cm^2/s] ~ 14.0e-5 / (0.89^1.1 * V_A^0.6) ~ 1.59e-4 / V_A^0.6
```

So **MW exponent ~ -0.6** — same family as Wilke-Chang and Hayduk-Laudie.
Numerical predictions agree with Wilke-Chang to within ~5–10% across
the typical drug-like MW range.

## Training set
Small aqueous compendium (~40 compounds) from Stokes, Hartley, and
related diffusion-cell experimental work in the 1940s.

## Reported performance
Original report claims agreement within ~5% over the training data;
modern re-evaluations (Hayduk & Laudie 1974) show ~8–15% mean error
on independent test sets.

## Validity / limitations
- Water-as-solvent only.
- Does not account for solvent association (no phi term as in W-C).
- LeBas-volume convention.
- Mostly displaced by Wilke-Chang (1955) and Hayduk-Laudie (1974)
  for current use.

## Notes
Included for historical completeness — the V_A^(-0.6) and eta^(-1.1)
exponent pair is essentially the same as Hayduk-Laudie's later refit.
Useful as a quick paper-back-of-envelope alternative when LeBas
volumes are pre-tabulated.

## References
- DOI 10.1021/ie50519a036 (Othmer & Thakar 1953).
- Hayduk W, Laudie H. AIChE J. 1974;20(3):611–615 (modern refit).
