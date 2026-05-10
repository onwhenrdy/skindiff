---
name: chen_lian_han_2010
predicts: D
layer: SC
units: cm^2/s
descriptors: MW (Da), log K_ow, h-bond donor/acceptor count
domain: hydrophilic + hydrophobic small solutes; bricks-and-mortar SC
status: modern
---

# Chen–Lian–Han 2010 — Bricks-and-mortar SC permeation with structural D

## Citation
Chen L, Lian G, Han L. "Modeling transdermal permeation. Part I. Predicting
skin permeability of both hydrophobic and hydrophilic solutes." AIChE J.
2010;56(5):1136–1146. DOI: 10.1002/aic.12048.
(Part II / extensions: AIChE J. 2012, J Pharm Sci 2013.)

## Equation
A two-dimensional bricks-and-mortar SC with three diffusivities:

- **Lipid lateral D_lip** — uses Mitragotri 2003 free-volume scaled-particle
  form:

```
D_lip(r) = D_lip_0 * exp(-pi * (r / r_f)^2)        (r in A, r_f ~ 4 A)
```

producing D_lip in the 6e-10 to 6e-7 cm^2/s range for MW 100–500 Da
(reported by Chen et al.).

- **Corneocyte D_cor** — Wilke-Chang baseline times Bungay-Brenner
  hindrance polynomial in lambda = r/r_fiber:

```
D_cor = D_aq * (1 - phi_f) * H(lambda)        H = (0.9999 - 1.276 lam + ...)
```

(same as `nitsche_kasting_2013_corneocyte`).

- **Aqueous-pore D_pore** — D_aq scaled by tortuosity factor 1/T (T ~ 2–4)
  for hydrophilic solutes through transient water channels.

The headline predictive result combines them: the **MW exponent for
log P_SC** comes out to slope ~ −0.0085 in the Chen et al. simple QSAR
overlay (logKp = −2.55 + 0.65 log P − 0.0085 MW, R^2 = 0.91 on their
test set) — so for log P_SC ~ log D + log K, the lipid D contributes
essentially the Mitragotri exponential and the corneocyte D the
Wilke-Chang power-law.

## Training set
Validated against 124 compounds from the Lian–Chen–Han 2008 curated
permeability database (also used by Lian 2008 model bake-off).
Roughly half hydrophobic (log K_ow > 1), half hydrophilic.

## Reported performance
RMSE log K_p ~ 0.6 across the database, comparable to or slightly
better than Mitragotri 2003 alone; clearly superior to Potts–Guy
for log K_ow < 0 hydrophilic solutes (where the aqueous-pore term
matters). R^2 ~ 0.85–0.91 depending on subset.

## Validity / limitations
- Inherits limitations of each component model.
- Predictions for ionised species are weak (no explicit charge handling).
- Geometric parameters (r_f, phi_f, tortuosity) are fixed by SC
  microstructure — no patient-to-patient variation built in.

## Notes
This paper is the cleanest reference for "use Mitragotri D in lipid +
Bungay-Brenner D in corneocyte" as a turn-key recipe. It directly
inspires the assembly used in Hansen et al. and in modern in silico
skin tools.

## References
- DOI 10.1002/aic.12048 (Chen, Lian, Han 2010 Part I).
- Lian G, Chen L, Han L. J Pharm Sci. 2008;97(1):584–598.
  DOI 10.1002/jps.21074 (model evaluation across 7 QSARs).
- Hansen S, Lehr C-M, Schaefer UF. Adv Drug Deliv Rev. 2013;65(2):251–264
  (assembly review).
