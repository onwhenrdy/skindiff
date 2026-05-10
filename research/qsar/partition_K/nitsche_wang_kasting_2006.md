---
name: nitsche_wang_kasting_2006
predicts: K
phase_pair: SC/water (whole tissue, additive two-phase decomposition)
units: dimensionless (volume-based; intracellular water set as reference)
descriptors: logKow
domain: 72 SC/water partition data points spanning more than 8 decades of log Kow (approx -3 to +6), MW ~ 18-500 Da
status: benchmark
---

# Nitsche, Wang & Kasting (2006) two-phase SC/water partition

## Citation
Nitsche JM, Wang TF, Kasting GB. A two-phase analysis of solute partitioning into the stratum corneum. J Pharm Sci. 2006 Mar;95(3):649-666. DOI: 10.1002/jps.20549. PMID: 16432875.

## Equation
Whole-tissue K_SC/w is decomposed by volume-fraction sum:

    K_SC/w = phi_lip * K_lip/w + (1 - phi_lip) * K_cor/w
    K_cor/w = phi_w_cor + phi_pro * PC_pro/w

with the best-fit power laws against log Kow:

    K_lip/w  = 0.43 * K_ow^0.81
    PC_pro/w = 4.2  * K_ow^0.31

For fully hydrated SC the volume fractions are approximately phi_lip ~ 0.13, phi_pro ~ 0.40, phi_w_cor ~ 0.47.

**K-vs-Kow exponent: lipid phase 0.81, protein phase 0.31, whole-tissue effective slope 0.5-0.8 depending on log Kow regime.** This is the most-cited recalibration of Anderson-Raykar's two-domain decomposition.

## Training set
72 K_SC/w experimental data points from the literature, spanning ~8 decades of K_ow. Compounds are mostly small organic permeants from Flynn-style permeability databases; subset that have an *independently measured* SC/water partition coefficient (not back-calculated from K_p).

## Reported performance
RMS error in log10 K_SC/w of 0.30, attributed by the authors to be at the floor of experimental scatter (i.e. essentially limited by the noise in the measurements).

## Validity / limitations
- Hydrated human SC.
- Designed and validated as the *whole-tissue* K_SC/w predictor through the additive sum -- the individual K_lip/w and PC_pro/w fits should not be applied piecewise unless the user separately models the lipid and corneocyte volume fractions correctly.
- Larger, more flexible molecules (peptides, MW > 500) are not in the training set.
- Ionizable compounds: use the neutral-form Kow at the relevant pH; the model does not handle ionization explicitly.

## Notes
For `skindiff`, this is the recommended source for SC-layer K when modeling SC as a single homogeneous slab. If SC is split into a lipid sub-layer and a corneocyte sub-layer (as some advanced uses of `skindiff` would), the two power laws apply to those sub-layers individually. The 0.81 exponent on the lipid phase is the most influential number in the literature for highly lipophilic compounds.

The 2006 paper supersedes Anderson-Raykar 1988-89 as the calibration of the two-domain framework; Anderson-Raykar provided the *form*, Nitsche-Wang-Kasting provided the *numbers* on the much larger pooled dataset.

## References
- Nitsche, Wang, Kasting 2006 doi:10.1002/jps.20549
- Source decomposition: Raykar, Fung, Anderson 1988 doi:10.1023/A:1015956705293
