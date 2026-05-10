---
name: hayduk_laudie_1974
predicts: D
layer: aqueous
units: cm^2/s
descriptors: solute molar volume V_A (LeBas), water viscosity eta_w
domain: dilute aqueous solutions of non-electrolytes; 10–30 C
status: classical
---

# Hayduk–Laudie 1974 — Aqueous diffusion (water-only specialization)

## Citation
Hayduk W, Laudie H. "Prediction of diffusion coefficients for nonelectrolytes
in dilute aqueous solutions." AIChE J. 1974;20(3):611–615.
DOI: 10.1002/aic.690200329.

## Equation
A water-specific simplification of Wilke–Chang, optimized on aqueous data:

```
D_AB [cm^2/s] = 13.26e-5 * eta_w^(-1.14) * V_A^(-0.589)
```

with:
- eta_w = water dynamic viscosity in centipoise (~0.89 cP at 25 C)
- V_A = solute molar volume in cm^3/mol (LeBas group additivity at
  normal boiling point)

Dropping the explicit T and using 25 C water (eta_w = 0.89 cP):

```
D_aq [cm^2/s] ~ 13.26e-5 * 0.89^(-1.14) * V_A^(-0.589)
            ~ 1.51e-4 / V_A^0.589
```

Since V_A ~ MW for organic compounds (rough heuristic), **effective MW
exponent ~ -0.589**, very close to Wilke–Chang's −0.6 (these two forms
differ by < 5% in D_aq predictions for most small organics).

For MW = 100, D_aq ~ 1.0e-5; MW = 500, D_aq ~ 4e-6 cm^2/s.

## Training set
87 organic compounds in water, MW 30–400, mostly alkanes, alcohols,
amino acids, sugars, drug-like molecules. Includes both polar and
non-polar solutes.

## Reported performance
Mean absolute error ~5.8% on the training set; ~10% on independent test
compounds. Slightly **better than Wilke–Chang for water as solvent** —
the dropped sqrt(phi*M_B) term is essentially constant for water and
absorbed into the 13.26e-5 prefactor.

## Validity / limitations
- Water as solvent only — not extensible to organic solvents.
- Non-electrolytes; ionized species need a Nernst-Einstein correction
  using ionic charge and mobility.
- LeBas additive volumes may need adjustment for highly branched or
  cyclic structures.
- Temperature range typically 10–40 C; the eta_w^(-1.14) explicit
  dependence handles this within range.

## Notes
This is the most defensible aqueous-D baseline for skin-diffusion
modeling at body T (37 C), since it avoids the awkward water
association factor phi = 2.6 of Wilke–Chang. Many recent
finite-element / multiphase skin codes (Hansen et al. assembly,
Frasch & Barbero) default to Hayduk–Laudie over Wilke–Chang.

## References
- DOI 10.1002/aic.690200329 (Hayduk & Laudie 1974).
- Tucker WA, Nelken LH. In *Handbook of Chemical Property Estimation
  Methods*, eds. Lyman, Reehl, Rosenblatt. McGraw-Hill, 1982 (often-
  cited summary tabulation).
