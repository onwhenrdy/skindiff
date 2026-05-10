---
name: lian_chen_han_2008_eval
predicts: D (indirect, via P)
layer: SC
units: cm/s for P; cm^2/s for inferred D
descriptors: MW, log K_ow (and model-specific descriptors)
domain: comparative benchmark across 7 SC permeability/D models on 124 compounds
status: benchmark
---

# Lian–Chen–Han 2008 — Evaluation of seven skin-permeability models

## Citation
Lian G, Chen L, Han L. "An evaluation of mathematical models for predicting
skin permeability." J Pharm Sci. 2008;97(1):584–598.
DOI: 10.1002/jps.21074. PMID 17696163.

## Equation
This is **not** itself a D QSAR but a head-to-head comparison of seven
skin-permeability / diffusion models against a curated benchmark dataset.
Key models compared:

1. Potts–Guy 1992 (log Kp = a + b log K_ow + c MW; not a D model)
2. Robinson–Riviere QSAR
3. Lien–Gao QSAR
4. Cleek–Bunge stratum corneum + viable epidermis two-layer model
5. **Mitragotri 2002** scaled-particle theory (D_lip ~ exp(-0.46 r^2))
6. Wang–Kasting–Nitsche 2007 multiphase microscopic (D_lip-lat ~ MW^-2.43)
7. Their own Chen–Lian–Han 2010 precursor

Inferred D values vary across models by 1–4 log decades for any given
solute, depending mostly on the MW exponent assumed.

## Training set / benchmark
124 compounds from in-vitro human skin permeation, MW 18–765 Da, log K_ow
range −2.7 to 5.5. Compiled from the Wilschut, Cronin, Vecchia datasets
and refined by the authors.

## Reported performance
Best mechanistic model: **Mitragotri 2002** (lowest MSE across the 124
compounds, log K_p RMSE ~ 0.55).
Best simple empirical model: Potts-Guy (RMSE ~ 0.65).
Wang-Kasting matches Mitragotri on aggregate but is more accurate on
hydrophilic solutes.

## Validity / limitations
- Not a D model on its own — converts to D via D = P*h/K with assumed h
  and K_lip.
- The benchmark dataset is dominated by lipophilic solutes (~80%),
  which biases the verdict.
- Re-evaluation on later, larger datasets (Magnusson 2004, Vecchia 2003
  consolidated) gives smaller margins between Mitragotri, Wang–Kasting,
  and Chen.

## Notes
Often cited as the canonical "Mitragotri wins" result. The relative
ranking between Mitragotri's exp(−c r^2) and Wang–Kasting's MW^(−2.43)
is sensitive to the MW > 300 tail of the dataset — both models are
plausible on the 100–300 Da bulk.

## References
- DOI 10.1002/jps.21074 (Lian, Chen, Han 2008).
- DOI 10.1002/jps.10031 (Mitragotri 2002, the winner).
- DOI 10.1002/jps.20883 (Wang–Kasting–Nitsche 2007).
