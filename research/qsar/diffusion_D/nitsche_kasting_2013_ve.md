---
name: nitsche_kasting_2013_ve
predicts: D
layer: VE
units: cm^2/s
descriptors: MW (via D_aq), molar volume V_A
domain: viable epidermis as cytoplasm + extracellular tortuous network
status: modern
---

# Nitsche–Kasting 2013 — Viable epidermis effective diffusion

## Citation
Nitsche JM, Kasting GB. "A microscopic multiphase diffusion model of viable
epidermis permeability." Biophys J. 2013;104(10):2307–2320.
DOI: 10.1016/j.bpj.2013.04.009. PMID 23708373; PMC3660641.

## Equation
The viable epidermis (VE) is modeled as a brick layer of keratinocytes
(cytoplasm) bathed in a thin extracellular fluid layer; both phases are
penetrable but each hinders solute relative to bulk water. The component
diffusivities are simple constant hindrance factors on the aqueous D:

```
D_cyt [cm^2/s] = D_aq / H_cyt          H_cyt = 3
D_ext [cm^2/s] = D_aq / H_ext          H_ext = 2
```

Tissue-scale homogenization gives an effective transverse diffusivity

```
D_VE_eff [cm^2/s] = D_aq / H_epi       H_epi = 5..10  (typical 7)
```

D_aq is the bulk water diffusivity of the solute, computed from
Wilke–Chang (`wilke_chang_1955`) or Hayduk–Laudie (`hayduk_laudie_1974`).
Thus the VE diffusivity scales with MW only via D_aq:

```
D_VE_eff ~ D_aq ~ V_A^(-0.6) ~ MW^(-0.6)         (loosely, since V_A ~ MW)
```

## Training set
Hindrance factors H_cyt, H_ext, H_epi calibrated against a small set
of canonical compounds: water, ethanol, nicotinamide, testosterone
(plus literature D_VE for hydrocortisone). No formal training-set
regression — these are "provisional estimates" from contemporary
knowledge per the authors.

## Reported performance
For water and hydrocortisone, model D_VE matches measured D within
factor ~2; for highly bound solutes (e.g. testosterone) the model
under-predicts the apparent D unless an effective binding factor is
added (cf. Kretsos–Kasting 2008 dermis model).

## Validity / limitations
- Hindrance factors are constant (size-independent) — the model
  ignores molecular-size-dependent steric hindrance in the cytoplasm,
  which begins to matter for solutes > ~500 Da.
- Does not separate intracellular vs. paracellular (tight-junction)
  routes, which can dominate for hydrophilic large solutes.
- All non-water-like behavior (binding, charge) is offloaded to K and
  to D_aq; D_VE itself has no fitted MW dependence.

## Notes
The simplicity is intentional: VE is not the rate-limiting layer for
most topical formulations, so a one- or two-parameter hindrance model
is adequate. For PBPK-type applications an effective layer thickness
~50–80 um with H_epi = 7 reproduces published K_p_VE within an order
of magnitude across the standard QSAR test sets.

## References
- DOI 10.1016/j.bpj.2013.04.009 (Nitsche & Kasting 2013).
- See also Kretsos K, Kasting GB. Skin Pharmacol Physiol 2005;18(2):55–74.
  DOI 10.1159/000083706 (capillary clearance, related VE/dermis modeling).
