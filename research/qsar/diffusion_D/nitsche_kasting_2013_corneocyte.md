---
name: nitsche_kasting_2013_corneocyte
predicts: D
layer: corneocyte
units: cm^2/s
descriptors: MW (via D_aq), molecular_radius_r, fiber-volume-fraction phi_f
domain: corneocyte (hydrated keratin matrix) interior of SC; small molecules
status: modern
---

# Nitsche–Kasting 2013 — Hindered diffusion in corneocyte / keratin matrix

## Citation
Nitsche JM, Kasting GB. "A microscopic multiphase diffusion model of viable
epidermis permeability." Biophys J. 2013;104(10):2307–2320.
DOI: 10.1016/j.bpj.2013.04.009. PMID: 23708373; PMC3660641.
(Companion: Hansen et al. on corneocyte transport, J Pharm Sci 2011/2013.)

## Equation
Corneocyte diffusion is computed from aqueous diffusion D_aq (Wilke-Chang
or Hayduk-Laudie) hindered by the keratin fiber matrix, using the
Brenner–Bungay-style fiber-hindrance polynomial:

```
D_cor = D_aq * (1 - phi_f) * (0.9999 - 1.2762*lambda + 0.0718*lambda^2 + 0.1195*lambda^3)
```

where:
- D_aq is bulk aqueous diffusivity (Wilke-Chang form, see `wilke_chang_1955`)
- phi_f is the keratin volume fraction (~0.55 for fully hydrated SC)
- lambda = r_solute / r_fiber, the ratio of solute radius to keratin
  fiber radius (r_fiber ~ 35–40 A)

For viable epidermis (no keratin), the same paper uses simpler hindrance
factors:

```
D_cyt = D_aq / H_cyt        H_cyt = 3       (cytoplasm)
D_ext = D_aq / H_ext        H_ext = 2       (extracellular fluid)
D_epi_eff ~ D_aq / H_epi    H_epi = 5..10   (whole-tissue effective)
```

The MW dependence enters only via D_aq (~MW^-0.5 to ~MW^-0.6 depending
on aqueous correlation). For water through fully-hydrated SC corneocytes,
D_cor / D_aq ~ 0.10–0.30, dropping to ~0.01–0.05 as solute approaches
fiber size.

## Training set
Provisional parameter set from "contemporary knowledge" (the authors'
phrasing) rather than a regressed fit. Validation against a small
"limited anecdotal" dataset of 4 representative compounds: water,
ethanol, nicotinamide, testosterone. Hansen et al.'s extended keratin
binding database (~100 compounds) backs the K_cor side of the model.

## Reported performance
The model gives D_epi_eff ~ 0.10–0.20 * D_aq, in agreement with
permeability of viable epidermis to water and hydrocortisone within
order of magnitude.

## Validity / limitations
- Hindrance polynomial assumes lambda < 0.95 (solute notably smaller
  than fiber). For macromolecules lambda → 1 the polynomial blows up.
- phi_f is hydration-dependent — in dry SC corneocytes (phi_f → 0.7),
  D_cor falls to near zero by the polynomial but biologically there
  is no transport at all. Use the model only at physiological hydration.
- Does not include solute-keratin binding (which slows transport
  beyond steric hindrance). For lipophilic solutes with high keratin
  affinity see Nitsche–Frasch 2022 framework.

## Notes
This is **the** standard model for corneocyte / viable-epidermis
diffusion in the modern multiphase SC literature (Wang–Kasting,
Hansen, Frasch). It explicitly separates the aqueous-baseline D_aq
(via Wilke–Chang) from the geometric hindrance — so the MW scaling
inherits from the aqueous correlation, not from any SC-specific fit.

## References
- DOI 10.1016/j.bpj.2013.04.009 (Nitsche & Kasting 2013).
- Hansen S, Naegel A, Heisig M et al. J Pharm Sci. 2011;100(5):1655–1668.
  DOI 10.1002/jps.22397 (extended keratin binding dataset).
- Bungay PM, Brenner H. Int J Multiphase Flow 1973 (hindrance polynomial).
