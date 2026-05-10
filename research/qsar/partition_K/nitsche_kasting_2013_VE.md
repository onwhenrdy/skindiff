---
name: nitsche_kasting_2013_VE
predicts: K
phase_pair: VE/water (viable epidermis vs water)
units: dimensionless (volume-based; intracellular water reference)
descriptors: logKow (lipid sub-fraction); MW for diffusion side
domain: small organic permeants used to validate the multiphase microscopic model; same domain as Wang-Kasting-Nitsche 2007 SC model
status: modern
---

# Nitsche & Kasting (2013) -- viable epidermis K_VE/water

## Citation
Nitsche JM, Kasting GB. A microscopic multiphase diffusion model of viable epidermis permeability. Biophys J. 2013 May 21;104(10):2307-2320. DOI: 10.1016/j.bpj.2013.03.056. PMC: PMC3660641.

## Equation
Volume-fraction-weighted average over cytoplasm, lipid, and extracellular fluid:

    K_VE/w  =  phi_cyt * K_cyt/w  +  phi_lip * K_lip/w  +  phi_ext * K_ext/w

with the standard hydrated-VE composition:

    phi_cyt ~ 0.879   (cytoplasm)
    phi_lip ~ 0.001   (lipid bilayer + organelle membrane)
    phi_ext ~ 0.120   (extracellular fluid)

K_cyt/w follows a corneocyte-like protein-binding correlation (close to the Anderson-Raykar / Hansen 0.31-0.32 keratin slope but with cytoplasmic protein density), K_lip/w follows the lipid bilayer correlation (Klip/w = 0.43 * Kow^0.81 from Wang/Kasting), and K_ext/w ~ 1 (interstitial fluid is essentially bulk-aqueous-like).

For most small solutes the lipid contribution phi_lip * K_lip/w is negligible because phi_lip is so small; the result is

    K_VE/w  approx  1   (within a factor of 2 for log Kow up to about +4)

For very lipophilic compounds (log Kow > 5) the tiny lipid fraction starts to lift K_VE/w above 1. **K-vs-Kow exponent for VE is essentially 0 across most of the practical range and tends toward the lipid exponent (~0.8) only at very high log Kow.**

## Training set
Calibrated using the same K_lip/w / K_cor/w correlations as Wang-Kasting-Nitsche 2007; volume fractions taken from Wang-Kasting cytoplasm composition data.

## Reported performance
The 2013 paper validates the model architecture against measured viable-epidermis permeabilities for a small panel of compounds. K-side performance is inherited from the underlying lipid and protein-binding correlations.

## Validity / limitations
- Hydrated viable epidermis.
- Cytoplasmic K_cyt/w uses an effective protein-binding constant; for compounds with specific cytosolic-protein binding (heat-shock proteins, drug-transporter substrates) the model under-predicts K_VE/w.
- Charged compounds need ionization correction at intracellular pH ~ 7.0-7.2.

## Notes
For `skindiff`, viable-epidermis K is typically of order 1 (similar to dermis), with a weak Kow dependence; this is the appropriate prior. If `skindiff` includes a viable-epidermis layer between SC and dermis, a flat K_VE = 1 is a defensible default for non-lipophilic compounds, and the Nitsche-Kasting 2013 volume-fraction sum is the right way to add a logKow correction for lipophilic ones.

## References
- Nitsche & Kasting 2013 doi:10.1016/j.bpj.2013.03.056
- K_lip/w source: Nitsche, Wang, Kasting 2006 doi:10.1002/jps.20549
- K_cyt/w source: Hansen et al. 2011 doi:10.1002/jps.22396
