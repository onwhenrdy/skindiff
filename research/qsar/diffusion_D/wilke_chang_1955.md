---
name: wilke_chang_1955
predicts: D
layer: aqueous
units: cm^2/s
descriptors: MW (solvent), molar volume V_A (solute, at boiling point), T, eta
domain: dilute aqueous (or organic) solution; non-electrolytes; ~10% accuracy
status: classical | benchmark
---

# Wilke–Chang 1955 — Aqueous diffusivity baseline

## Citation
Wilke CR, Chang P. "Correlation of diffusion coefficients in dilute solutions."
AIChE J. 1955;1(2):264–270. DOI: 10.1002/aic.690010222.

## Equation
The most widely used aqueous-diffusivity correlation:

```
D_AB [cm^2/s] = 7.4e-8 * (phi * M_B)^0.5 * T / ( eta_B * V_A^0.6 )
```

with:
- phi = solvent association factor (water: 2.6 in original; 2.26 also seen)
- M_B = solvent molecular weight (water: 18.0 g/mol)
- T = temperature in Kelvin
- eta_B = solvent dynamic viscosity in centipoise (water at 25 C: 0.890 cP)
- V_A = solute molar volume at normal boiling point in cm^3/mol
  (LeBas additive estimate)

For aqueous diffusivity at 25 C this collapses to

```
D_aq [cm^2/s] ~ 7.4e-8 * sqrt(2.6 * 18) * 298 / (0.890 * V_A^0.6)
            ~ 1.95e-4 / V_A^0.6
```

Since V_A scales roughly linearly with MW for organic compounds (typical
density 1 g/mL → V_A ≈ MW), the **effective MW exponent is approximately
−0.6**:

```
D_aq ~ MW^(-0.6)        (approximation when V_A ~ MW)
```

Distinct from Stokes–Einstein D ~ MW^(-1/3) and from Hayduk–Laudie's
slightly steeper V_A^(-0.589). For MW = 100, D_aq ~ 1.2e-5 cm^2/s; for
MW = 500, D_aq ~ 5e-6 cm^2/s. A shift in scaling form by Hansen et al.
gives the alternative piecewise:

```
D_aq = 1.92e-4 / V_A^0.6        if V_A <= 445.2 cm^3/mol
D_aq = 3.78e-5 / V_A^(1/3)      if V_A >  445.2 cm^3/mol
```

(the high-V_A branch matches Stokes–Einstein for large solutes).

## Training set
~250 solute-solvent pairs across water and organic solvents, MW 18–400,
mostly hydrocarbons, alcohols, and small organics. Compiled from earlier
experimental compendia.

## Reported performance
Stated accuracy ±10% for non-associating, non-electrolyte solutes in
dilute solution. Outliers up to ±30% for hydrogen-bonded systems and
near critical conditions. Significantly better for water-as-solvent
than for organic-as-solvent.

## Validity / limitations
- Dilute solution only.
- Non-electrolytes only.
- Solute molar volume at the *boiling point* (LeBas group additivity);
  not the same as molar volume at room T.
- Association factor phi for water is empirical; original Wilke/Chang
  recommend 2.6 but 2.26 has been re-fit by some workers.
- Underestimates D for very small solutes (water in water: predicted
  2.0e-5, actual 2.3e-5 at 25 C).

## Notes
Used as the **D_aq baseline** by all SC and VE multiphase models
(Nitsche–Kasting, Wang–Kasting, Hansen) — D_cor and D_VE are computed
as hindrance factor times Wilke–Chang. So the −0.6 V_A exponent
propagates indirectly into all skin diffusion predictions, regardless
of which SC model is on top.

## References
- DOI 10.1002/aic.690010222 (Wilke & Chang 1955).
- Reid RC, Prausnitz JM, Poling BE. *The Properties of Gases and Liquids*,
  4th–5th ed., McGraw-Hill (chapter on liquid diffusion).
