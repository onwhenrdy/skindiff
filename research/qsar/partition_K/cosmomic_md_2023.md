---
name: cosmomic_md_2023
predicts: K
phase_pair: SC_lipid/water
units: dimensionless (volume-based)
descriptors: full molecular structure (no summary descriptor); COSMOmic uses sigma profiles, MD uses force-field-resolved atoms
status: modern (in-silico, not a closed-form QSAR)
---

# COSMOmic + MD predictions of K_lip/water (Piasentin 2023 and related)

## Citation
Piasentin N, Lian G, Cai Q. In silico prediction of stratum corneum partition coefficients via COSMOmic and molecular dynamics simulations. J Phys Chem B. 2023 Apr 6;127(13):2719-2729. DOI: 10.1021/acs.jpcb.2c08566. PMID: 36930176. PMC: PMC10068742.

## Equation
This is not a single closed-form QSAR. COSMOmic computes an inhomogeneous solubility profile across a model SC lipid bilayer using COSMO-RS sigma profiles and a precomputed "membrane sigma profile" for the SC bilayer, then integrates the Boltzmann factor across the bilayer to get K_lip/w. MD computes K_lip/w by free-energy umbrella sampling across a CHARMM/GROMACS atomistic SC bilayer model.

When the COSMOmic / MD predictions are summarised as a regression against log Kow on the test set used in this paper:

    log K_lip = alpha + beta * log K_ow,
    beta = 0.74, RMSE = 0.46, R = 0.84       (best-fit through the experimental data they used)

with cited literature values from earlier QSARs of beta = 0.69 (Frasch-Barbero / Wang), 0.70 (Mitragotri / Potts-Guy implicit), 0.86 (Johnson 1996). **K-vs-Kow exponent across the in-silico predictions is in the same 0.7-0.8 band as the empirical fits**, confirming that the lipid-bilayer / octanol slope in the SC is robustly less than 1 and centred near 0.7.

## Training set
Test set of 16-20 compounds with experimental K_lip/w from the Hansen 2013 curated dataset and earlier liposome partitioning measurements; not "training" in the traditional sense (the COSMOmic / MD predictions are de novo from molecular structure, then compared to experiment).

## Reported performance
RMSE in log10 K_lip/w of ~0.46 against experimental measurements; R = 0.84. Comparable to the empirical QSARs but with the advantage that no Kow descriptor is required as input -- only the molecular structure.

## Validity / limitations
- Computationally expensive (especially MD) -- not suitable for high-throughput screening.
- Force-field / sigma-profile choice matters; different SC bilayer models give noticeably different K predictions.
- Does not address keratin/protein binding; only the lipid-phase K.

## Notes
For `skindiff`, COSMOmic / MD K_lip/w predictions are useful when log Kow is unknown or unreliable (e.g. novel compounds, compounds where Kow is hard to measure due to surfactant behaviour). For routine use, the closed-form Hansen 2013 / Wang-Kasting 2007 K_lip/w correlations are faster and within their error bars consistent with the COSMOmic predictions.

## References
- Piasentin, Lian, Cai 2023 doi:10.1021/acs.jpcb.2c08566
- COSMOmic methodology: Klamt et al. earlier; refer to the Piasentin 2023 reference list.
