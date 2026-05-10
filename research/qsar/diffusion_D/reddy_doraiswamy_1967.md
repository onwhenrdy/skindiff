---
name: reddy_doraiswamy_1967
predicts: D
layer: aqueous
units: cm^2/s
descriptors: solute molar volume V_A, solvent molar volume V_B, T, eta
domain: dilute liquid solutions, water and organic solvents; non-electrolytes
status: classical
---

# Reddy–Doraiswamy 1967 — Generalized Stokes–Einstein-style aqueous D

## Citation
Reddy KA, Doraiswamy LK. "Estimating Liquid Diffusivity." Ind Eng Chem
Fundam. 1967;6(1):77–79. DOI: 10.1021/i160021a012.

## Equation
A semi-empirical extension of Stokes–Einstein with explicit solute and
solvent molar-volume dependence:

```
D_AB [cm^2/s] = K * (M_B^0.5) * T / [ eta_B * (V_A * V_B)^(1/3) ]
```

with the constant K piecewise on V_A / V_B ratio:
- K = 10e-8 if V_B / V_A < 1.5
- K = 8.5e-8 if V_B / V_A >= 1.5

units: M_B in g/mol, T in K, eta_B in cP, V_A and V_B in cm^3/mol.

For water (V_B = 18.07, M_B = 18, eta_B = 0.89 cP at 25 C):

```
D_aq [cm^2/s] ~ K * sqrt(18) * 298 / (0.89 * (V_A * 18.07)^(1/3))
            ~ 1.42e-3 / V_A^(1/3)        (K=10e-8 small solutes)
```

So **effective MW exponent = -1/3**, matching Stokes–Einstein. This is
shallower than Wilke–Chang's −0.6; consequently Reddy–Doraiswamy
*overpredicts* D for large solutes relative to W–C and Hayduk–Laudie.

## Training set
~96 binary solute-solvent systems compiled from earlier diffusivity
literature. Mix of aqueous and organic solvents; small organics.

## Reported performance
Average percent error ~13–18% on the original training set; comparable
to but slightly worse than Wilke–Chang for water, slightly better for
some organic solvents.

## Validity / limitations
- Non-electrolytes only.
- The piecewise K constant is awkward and creates a small discontinuity
  at V_B/V_A = 1.5.
- LeBas molar volumes assumed.
- Domain weakly defined for solutes >300 cm^3/mol.

## Notes
Less commonly used than Wilke–Chang or Hayduk–Laudie in current skin
modeling; included here mostly to document the **MW^(-1/3) family**
of aqueous correlations against the **MW^(-0.6)** Wilke-Chang line.
The two predictions can differ by 2–3x at MW ~ 500.

## References
- DOI 10.1021/i160021a012 (Reddy & Doraiswamy 1967).
- Reid RC, Prausnitz JM, Poling BE. *The Properties of Gases and
  Liquids*, 5th ed., chapter 11.
