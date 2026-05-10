---
name: wang_kasting_2007_ktrans
predicts: D
layer: SC_lipid
units: cm/s (mass-transfer coefficient, not D)
descriptors: MW (Da)
domain: trans-bilayer hopping in SC; macroscopic-effective resistance
status: modern
---

# Wang–Kasting–Nitsche 2007 — Trans-bilayer mass-transfer coefficient

## Citation
Wang T-F, Kasting GB, Nitsche JM. "A multiphase microscopic diffusion model
for stratum corneum permeability. II. Estimation of physicochemical
parameters, and application to a large permeability database."
J Pharm Sci. 2007;96(11):3024–3051. DOI 10.1002/jps.20883.

## Equation
Bilayer-to-bilayer hopping is treated as a single-step mass transfer with
rate constant k_trans (units cm/s — a velocity), **not** a continuous D.
Wang fitted log k_trans to the SC permeability database with a Kasting-style
free-area form linear in cube-root MW:

```
log10 k_trans [cm/s] = -0.725 - 0.792 * MW^(1/3)            (MW in Da)
```

For MW = 100: k_trans ~ 5.6e-5 cm/s. For MW = 500: k_trans ~ 1.3e-7 cm/s.

The relationship to a continuous D is k_trans ~ D_trans / d_lamella, where
d_lamella is the inter-bilayer spacing (~5 nm). Thus implicitly

```
D_trans [cm^2/s] ~ d_lamella * k_trans = 5e-7 * 10^(-0.725 - 0.792 MW^(1/3))
```

The MW^(1/3) form is identical in spirit to Stokes–Einstein (D ~ 1/r ~
MW^(-1/3)) but in **log-linear** form rather than power-law: the activation
*energy* for crossing a bilayer scales linearly with the molecular linear
size.

## Training set
The k_trans regression is the *one* free parameter of the Wang–Kasting
2007 multiphase fit. Trained against ~120 compounds (Flynn / Vecchia /
Magnusson curated SC permeability data) holding the other four transport
parameters (D_cor, K_cor/w, D_lip-lat, K_lip/w) at independently-derived
values.

## Reported performance
Within the multiphase model the residual log K_p RMS ~0.5–0.6. Removing
or perturbing the k_trans formula by ±0.5 in the prefactor changes
predicted log Kp by similar magnitude.

## Validity / limitations
- This is **not** a diffusion coefficient but a hopping mass-transfer
  coefficient. The conversion to D_trans needs an assumed lamellar
  spacing.
- k_trans is the only fit parameter in the Wang–Kasting model; all
  the apparent extrapolative power lives in this regression.
- Domain: organic solutes 100–500 Da; not validated for water (MW = 18,
  H-bonding outlier) or large drugs (MW > 600).

## Notes
The MW^(1/3) form ties Wang–Kasting back to Kasting 1992 — Kasting's
log10 D = a − b*MW reduces to log10 D = a − b*MW^1 instead of MW^(1/3),
which is one of the few empirically distinguishable predictions between
the two SC models.

## References
- DOI 10.1002/jps.20883 (Wang–Kasting–Nitsche 2007, Part II).
- DOI 10.1002/jps.20509 (Wang–Kasting–Nitsche 2006, Part I).
- Kasting GB, Smith RL, Anderson BD. *Prodrugs* (1992) — original
  MW-linear-in-log form.
