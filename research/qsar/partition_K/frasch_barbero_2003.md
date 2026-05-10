---
name: frasch_barbero_2003
predicts: K
phase_pair: SC_lipid/water and SC_protein/water (decomposed)
units: dimensionless (volume-based, water reference)
descriptors: logKow
domain: small-organic permeants used in finite-element SC transport simulations; calibrated on Anderson-Raykar / Wang dataset; log Kow approx -2 to +5
status: classical (used in many subsequent FEM SC models)
---

# Frasch & Barbero (2003) -- two-phase SC partitioning for FEM models

## Citation
Frasch HF, Barbero AM. Steady-state flux and lag time in the stratum corneum lipid pathway: results from finite element models. J Pharm Sci. 2003 Nov;92(11):2196-2207. DOI: 10.1002/jps.10466. PMID: 14603505.

## Equation
The K parameterization used as input to the brick-and-mortar finite element transport models is the two-power-law decomposition

    PC_lip/w = K_ow^0.69
    PC_pro/w = 4.2 * K_ow^0.31

so the **K-vs-Kow exponent is 0.69 for the lipid phase and 0.31 for the protein/keratin phase**. These two relations are close but not identical to the Nitsche-Wang-Kasting (2006) recalibrations (0.81 / 0.31 with prefactor 0.43 on the lipid phase): Frasch-Barbero used a slightly older calibration of the lipid slope, with prefactor essentially 1 (no decade-scale shift between lipid and octanol).

## Training set
Inherited K parameterization from Anderson-Raykar 1988 / 89 dataset and from Mitragotri's lipid bilayer partitioning estimates; not re-fit from scratch in the 2003 paper. The 2003 paper's contribution is the FEM transport calculation, with the K-fit treated as input.

## Reported performance
The K-vs-Kow correlations were not refit -- performance on K_SC/w is inherited from the original Anderson-Raykar et al. fits (rms log10 K_SC/w ~ 0.4). The 2003 paper reports the FEM result that the steady-state flux through brick-and-mortar SC, with these K values, matches one-dimensional analytical predictions to within ~10% for typical aspect ratios.

## Validity / limitations
- Same restrictions as the Anderson-Raykar source: hydrated SC, log Kow [-2, +5], MW < 500 Da.
- The 0.69 vs Wang/Kasting 0.81 difference matters: at log Kow = +5 the two lipid-phase predictions differ by a factor of ~3-4, which is the main source of disagreement between the Frasch-Barbero and Wang-Kasting parameterizations downstream.

## Notes
The Frasch-Barbero K parameterization is the most-used one in *transport-side* SC modelling papers from approximately 2003 to 2010 (papers that use FEM brick-and-mortar geometries), while Nitsche-Wang-Kasting 2006/2007 dominates the *partition-side* fitting literature from 2006 onward. For `skindiff`, both are reasonable choices; Wang/Kasting is the more recent and larger-dataset calibration of the same model form.

## References
- Frasch & Barbero 2003 doi:10.1002/jps.10466
- Source K-fit: Anderson, Higuchi, Raykar 1988 doi:10.1023/A:1015989929342
