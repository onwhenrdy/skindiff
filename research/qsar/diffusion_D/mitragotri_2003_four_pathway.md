---
name: mitragotri_2003_four_pathway
predicts: D
layer: SC_lipid
units: cm^2/s
descriptors: molecular_radius_r (Angstrom), K_ow, MW
domain: hydrophilic + hydrophobic small solutes; SC, four parallel pathways
status: classical
---

# Mitragotri 2003 — Four-pathway SC permeation (free-volume D in lipid)

## Citation
Mitragotri S. "Modeling skin permeability to hydrophilic and hydrophobic
solutes based on four permeation pathways." J Control Release. 2003
Feb 14;86(1):69–92. DOI: 10.1016/S0168-3659(02)00321-8.

## Equation
P_total = P_lateral + P_pore + P_shunt + P_freeVolume.
The free-volume diffusion through lipid bilayers uses scaled-particle theory
with bilayer free volume of typical radius r_f ~ 4 A:

```
D_freeVolume(r) ~ D_0 * exp(-pi * (r / r_f)^2)         (r, r_f in Angstrom)
```

This is the same exponential-in-area form as the 2002 paper (the SPT
prefactor on r^2 in the exponent is pi/r_f^2 ~ 0.20 with r_f = 4 A;
fitted prefactors in 0.2–0.5 range have been reported).

The lateral-diffusion contribution uses Saffman–Delbrück scaling
D_lateral ~ k_B T / (4 pi mu_membrane * h) * [ln(mu_m h / (mu_w r)) − gamma],
which is **logarithmic in r** rather than exponential — contrasting
with the size dependence of the other pathways.

## Training set
Re-fit on Flynn's human-skin permeability database (~120 compounds);
hydrophilic solutes (caffeine, mannitol, sucrose, urea) included to test
the aqueous-pore term that was absent in the 2002 hydrophobic-only model.

## Reported performance
RMS log10(P) error ~ 0.7 across the full hydrophilicity range; matches
hydrophilic-solute permeability that the 2002 single-pathway model
under-predicts by 1–3 log decades.

## Validity / limitations
The free-volume contribution dominates for medium-size hydrophobic
solutes (MW 100–300, log K_ow > 1). Pore and shunt terms have larger
parameter uncertainty and are most relevant for hydrophilic /
high-MW solutes. The sub-models share fitted constants — using one D
expression in isolation outside the integrated framework can over-
or under-shoot by an order of magnitude.

## Notes
This paper introduced the now-canonical decomposition of SC transport
into four parallel routes; later multiphase models (Wang–Kasting,
Nitsche–Kasting) reorganize the same physics into bricks-and-mortar
geometry rather than parallel resistances.

## References
- DOI 10.1016/S0168-3659(02)00321-8 (Mitragotri 2003).
- Mitragotri S. Pharm Res. 2001;18(8):1018–1023 (in-situ K, D in SC lipids).
