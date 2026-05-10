---
name: wang_kasting_nitsche_2007
predicts: K
phase_pair: SC_lipid/water and corneocyte/water (microscopic phases)
units: dimensionless (volume-based)
descriptors: logKow (primary), MW (for diffusion side; not the K side directly)
domain: ~150 compounds in the consolidated SC permeability database, log Kow approximately -3 to +6, MW < 500 Da
status: benchmark
---

# Wang, Kasting & Nitsche (2006/2007) -- multiphase microscopic SC model, K parameterization

## Citation
- Wang TF, Kasting GB, Nitsche JM. A multiphase microscopic diffusion model for stratum corneum permeability. I. Formulation, solution, and illustrative results for representative compounds. J Pharm Sci. 2006 Mar;95(3):620-648. DOI: 10.1002/jps.20509. PMID: 16447176.
- Wang TF, Kasting GB, Nitsche JM. A multiphase microscopic diffusion model for stratum corneum permeability. II. Estimation of physicochemical parameters, and application to a large permeability database. J Pharm Sci. 2007 Nov;96(11):3024-3051. DOI: 10.1002/jps.20883. PMID: 17876780.

## Equation
Microscopic-scale partition coefficients used in the brick-and-mortar SC model:

    K_lip/w = 0.43 * K_ow^0.81             (carried over from Nitsche-Wang-Kasting 2006)
    K_cor/w = phi_w_cor + phi_pro * 4.2 * K_ow^0.31

with corneocyte composition phi_w_cor ~ 0.55, phi_pro ~ 0.45 in the standard "fully hydrated" parameterization. The corneocyte K is therefore approximately

    K_cor/w  approx  0.55 + 1.9 * K_ow^0.31

i.e. **K_cor-vs-Kow exponent is 0.31 in the lipophilic limit and approaches 1 (a flat constant) in the hydrophilic limit** because the aqueous fraction dominates.

## Training set
Used the 72-point K_SC/w dataset of Nitsche-Wang-Kasting 2006 to set K_lip/w and PC_pro/w; then validated against ~150 compound permeability database (Flynn / Vecchia-Bunge / Magnusson) by predicting K_p without further fitting K -- only the trans-bilayer mass-transfer coefficient k_trans was fit to permeability data.

## Reported performance
On the permeability prediction (which depends on both K and D), the model "somewhat" outperforms Potts-Guy in overall RMSE and dramatically outperforms it on lag-time prediction. K-only performance is inherited from the 2006 two-phase paper (rms log10 K_SC/w error 0.30).

## Validity / limitations
- Same domain restrictions as Nitsche-Wang-Kasting 2006.
- Corneocyte composition is treated as a fixed bulk-tissue average; the model does not address depth-dependent SC hydration gradients.
- Designed for a brick-and-mortar transport simulator -- the K_lip/w and K_cor/w numbers are intended as inputs to a multiphase diffusion model, not as a single composite K_SC/w to be plugged into a one-slab approximation (use Nitsche-Wang-Kasting 2006 K_SC/w for that).

## Notes
For a `skindiff` user who wants to split the SC into a lipid sub-layer and a corneocyte sub-layer with separate K and D, this is the cleanest published parameterization. The exponents (0.81 for lipid, 0.31 for protein) come straight from the 2006 paper; the contribution of this 2007 paper is the demonstration that the *transcellular* (corneocyte-traversing) pathway, not the intercellular lipid-only pathway, dominates K_p for most compounds when these K values are used together with a proper k_trans.

## References
- Wang, Kasting, Nitsche 2006 (Part I, formulation) doi:10.1002/jps.20509
- Wang, Kasting, Nitsche 2007 (Part II, parameterization) doi:10.1002/jps.20883
- Nitsche, Wang, Kasting 2006 (K calibration) doi:10.1002/jps.20549
