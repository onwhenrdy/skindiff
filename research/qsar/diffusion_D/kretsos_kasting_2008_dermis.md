---
name: kretsos_kasting_2008_dermis
predicts: D
layer: dermis
units: cm^2/s
descriptors: MW^(1/3), binding factor B (from K_de or protein-binding fraction)
domain: mammalian dermis (human, guinea pig, rat, mouse); MW 18–476 Da
status: modern | benchmark
---

# Kretsos–Miller–Zamora-Estrada–Kasting 2008 — D in mammalian dermis

## Citation
Kretsos K, Miller MA, Zamora-Estrada G, Kasting GB. "Partitioning, diffusivity
and clearance of skin permeants in mammalian dermis." Int J Pharm.
2008;346(1–2):64–79. DOI: 10.1016/j.ijpharm.2007.06.020. PMID 17703903.

## Equation
Dermis diffusivity is modeled as a free aqueous diffusivity (size-dependent)
attenuated by a binding factor that captures slowdown by extravascular
serum proteins:

```
D_de = D_free(MW) / (1 + B)
```

with the free-water term fitted against a liquid-like cube-root-MW
correlation:

```
log10 D_free [cm^2/s] = c_0 - c_1 * MW^(1/3)         (MW in Da)
```

(c_1 captures the Stokes-Einstein-like size dependence, equivalent in
form to Wang–Kasting's k_trans and to free-area theories — see
`wang_kasting_2007_ktrans`). The binding factor is

```
B = (K_de - K_water) / K_water
```

i.e. the same partition-derived attenuation Kretsos uses for K. Equivalent
to D_de = D_aq * f_unbound where f_unbound is the unbound fraction in
dermis. For very low protein binding D_de → D_aq; for high binding
(testosterone etc.) D_de can be 5–20x lower than D_aq.

## Training set
26 compounds with literature in vitro permeation through mammalian dermis,
MW 18–476 Da, across 4 species (human, guinea pig, rat, mouse). Validated
with in-house tritiated testosterone partition+diffusion measurements.

## Reported performance
log10 D_de RMSE ~ 0.5–0.7; predicted-vs-observed regression slope ~0.9.
Reproduces dermal capillary clearance decay lengths of 210–920 um with
clearance rate constants 0.9e-4 to 14e-4 1/s (Kretsos & Kasting 2005,
DOI 10.1159/000083706).

## Validity / limitations
- "Liquid-like" MW^(1/3) form is empirical for the free term; the
  binding-factor structure is the more important novelty.
- Requires K_de or protein-binding fraction as auxiliary input — D_de
  predictions are coupled to K_de QSARs (e.g. ionization, log K_ow).
- Calibrated on solutes ≤ 500 Da; macromolecules and biologics are
  out of domain (different transport regime through interstitium).
- Does not explicitly model collagen-fiber tortuosity (lumped into
  the prefactor).

## Notes
This is the standard *D in dermis* baseline used by PBPK skin models
(MoBi, GastroPlus dermal). Often paired with the dermal capillary
clearance Pe (perfusion) of the same series (Kretsos & Kasting 2005).

## References
- DOI 10.1016/j.ijpharm.2007.06.020 (Kretsos et al. 2008).
- DOI 10.1159/000083706 (Kretsos & Kasting 2005, capillary clearance).
- Ibrahim R, Nitsche JM, Kasting GB. J Pharm Sci. 2010;99(12):4928–4939.
  DOI 10.1002/jps.22216 (improved measurement method, ~30 compounds).
