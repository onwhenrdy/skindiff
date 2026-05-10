---
name: kretsos_kasting_2008
predicts: K
phase_pair: dermis/water (also viable epidermis/water by extension)
units: dimensionless (volume-based, free aqueous concentration as reference)
descriptors: logKow, fraction unbound in plasma (fu), pKa (for ionizable solutes), MW
domain: 26 compounds, MW 18-476 Da, four species (human, guinea pig, rat, mouse); log Kow ~ -3 to +6
status: classical (the canonical dermis K model)
---

# Kretsos, Miller, Zamora-Estrada & Kasting (2008) -- dermis K model

## Citation
Kretsos K, Miller MA, Zamora-Estrada G, Kasting GB. Partitioning, diffusivity and clearance of skin permeants in mammalian dermis. Int J Pharm. 2008 Jan 4;346(1-2):64-79. DOI: 10.1016/j.ijpharm.2007.06.020. PMID: 17703903.

## Equation
The dermis/water partition coefficient is built up from three contributions: free aqueous concentration, binding to extravascular plasma-like proteins (predominantly serum albumin), and partitioning into a small lipid sub-fraction:

    K_de/w  =  (1 - phi_lip - phi_pro) * f_neutral / fu_intra
              + phi_pro * K_pro/w
              + phi_lip * K_lip/w

with phi_lip ~ 0.01-0.02 (dermis lipid volume fraction is small) and phi_pro the protein volume fraction (mostly extravascular albumin equivalent). Practically, for non-ionizable, non-lipophilic solutes the answer collapses to roughly

    K_de/w  approx  0.7 to 0.9    (close to unity; dermis behaves like a slightly protein-loaded aqueous gel)

with a logKow-driven correction that becomes important only for highly lipophilic compounds (log Kow > 3), where the small phi_lip * K_lip/w term begins to dominate. **The K-vs-Kow exponent for dermis is therefore close to 0 over most of the practical log Kow range and approaches the lipid-phase exponent (~0.7-0.8) only for highly lipophilic compounds.**

The albumin-binding piece uses fu (the fraction of solute unbound in dilute plasma) as a per-compound experimental input or as a separately predicted quantity (Lobell-Sivars / Watanabe correlations). For testosterone, fu = 0.25 +/- 0.02 in 2% HSA was their reference measurement.

## Training set
26 compounds, MW 18-476 Da, four species. In-vitro permeation data from the literature combined with their own (3)H-testosterone permeation experiments.

## Reported performance
The dermis K model fits the in-vitro flux data well across the 26 compounds; the dominant uncertainty is in fu (plasma-protein-binding fraction), not in the model form itself.

## Validity / limitations
- Hydrated dermis only (the model assumes free water is bulk-aqueous-like).
- Highly protein-bound compounds need an experimental fu input -- predicted fu from QSAR is the largest source of error.
- The lipid-phase contribution uses the same K_lip/w correlation as the SC models (Mitragotri / Wang / Kasting style); the small phi_lip in dermis means the lipid-phase exponent on Kow only matters for log Kow > 3.

## Notes
For `skindiff`, dermis K is typically close to 1 (not the dramatically large or small numbers seen in SC). When fitting `skindiff` to permeation/penetration data, prior K_dermis ~ 0.7 with a weak logKow correction is the right starting point and a good prior. The Kretsos 2008 model is the standard reference for "dermis K is 0.7-0.9 plus a small lipid bump for lipophilic compounds."

A companion model for viable-epidermis K_VE/w is given by Nitsche & Kasting (2013) "A microscopic multiphase diffusion model of viable epidermis permeability" (Biophys J 104:2307; doi:10.1016/j.bpj.2013.03.056), with K_VE/w as a volume-weighted sum over cytoplasm (~88%), lipid (~0.1%), and extracellular fluid (~12%). For most solutes K_VE/w is also of order 1.

## References
- Kretsos, Miller, Zamora-Estrada, Kasting 2008 doi:10.1016/j.ijpharm.2007.06.020
- Companion VE model: Nitsche & Kasting 2013 doi:10.1016/j.bpj.2013.03.056
