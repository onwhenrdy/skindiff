---
name: wang_kasting_2007_dlip
predicts: D
layer: SC_lipid
units: cm^2/s
descriptors: MW (Da)
domain: SC lipid lamellae, lateral diffusivity in continuous lipid phase
status: modern
---

# Wang–Kasting–Nitsche 2007 — Lateral D in SC lipid bilayers (multiphase, II)

## Citation
Wang T-F, Kasting GB, Nitsche JM. "A multiphase microscopic diffusion model
for stratum corneum permeability. II. Estimation of physicochemical
parameters, and application to a large permeability database."
J Pharm Sci. 2007;96(11):3024–3051. DOI: 10.1002/jps.20883. PMID 17876780.
(Part I: J Pharm Sci. 2006;95(3):620–648. DOI 10.1002/jps.20509.)

## Equation
Lateral diffusivity in the continuous lipid phase of SC, as fit by
Wang–Kasting and re-arranged by Hansen et al. (2013):

```
D_lip-lat [cm^2/s] = 8.98e-3 * MW^(-2.43) + 2.34e-9        (MW in Da)
```

i.e. a steep **MW^-2.43 power law** with an additive 2.3e-9 cm^2/s floor —
unlike the *exponential-in-area* form of Mitragotri 2002. For MW = 100,
this gives D_lip-lat ~ 2.5e-7 cm^2/s; for MW = 500, ~3e-9 cm^2/s. The floor
matters above MW ~ 400, where the power-law term dies away.

Trans-bilayer hopping is treated separately via the **mass transfer
coefficient k_trans** rather than a diffusivity — fit in Wang 2006 as

```
log10 k_trans [cm/s] = -0.725 - 0.792 * MW^(1/3)            (MW in Da)
```

This is a Kasting-style **MW^(1/3) free-area** form (cube root of MW =
proportional to the molecular linear dimension). One trans-bilayer hop
is treated as a discrete event, not a continuous D.

Corneocyte diffusion D_cor is computed from D_aq via hindered diffusion
theory (Renkin/Bungay-Brenner) — see `nitsche_kasting_2013_corneocyte`.

## Training set
Part II's 4-parameter model fit to ~120 compounds from a curated SC
permeability database; D_lip-lat itself was constrained from independent
SC-stripping and FRAP data, k_trans was the one fit-only parameter.

## Reported performance
log10 K_p AAD ~ 0.4–0.6 across the database; the model improves on
Potts–Guy by separating lipid and protein phases, but not dramatically
on aggregate Kp metrics.

## Validity / limitations
- Lateral D is for the continuous lipid mortar of bricks-and-mortar SC,
  *not* trans-bilayer transport (covered by k_trans).
- The MW exponent of -2.43 is much steeper than Mitragotri's effective
  area exponent — caused by including very small (water, MW=18) and
  hydrophilic permeants in the fit, which pull D up at low MW.
- Floor 2.3e-9 cm^2/s is *not* physical free-volume limit; it's an
  empirical asymptote.

## Notes
This is the lateral D used inside the "bricks-and-mortar" multiphase model;
when collapsed to a single effective D_SC it gives ratios of 1:30 to 1:1000
to Mitragotri 2002 depending on MW. Use both to bound the uncertainty.

## References
- DOI 10.1002/jps.20883 (Wang–Kasting–Nitsche 2007, Part II).
- DOI 10.1002/jps.20509 (Wang–Kasting–Nitsche 2006, Part I).
- Hansen S, Naegel A, Heisig M, Wittum G, Neumann D, Kostka K-H, Meingassner
  J, Schaefer UF, Lehr C-M. Skin Pharmacol Physiol 2008/2013 mathematical
  model assembly review.
