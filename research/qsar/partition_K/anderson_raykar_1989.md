---
name: anderson_raykar_1989
predicts: K
phase_pair: SC/water (decomposed into SC_lipid/water + corneocyte/water)
units: dimensionless (volume-based concentration ratio at equilibrium)
descriptors: logKow
domain: hydrocortisone esters and a series of hydrophilic-to-lipophilic permeants used to fit the two-domain model; log Kow approximately -1 to +5
status: classical
---

# Anderson, Raykar & coworkers (1988-1989) two-domain SC partition model

## Citation
- Raykar PV, Fung MC, Anderson BD. The role of protein and lipid domains in the uptake of solutes by human stratum corneum. Pharm Res. 1988 Mar;5(3):140-150. DOI: 10.1023/A:1015956705293. PMID: 3244625.
- Anderson BD, Higuchi WI, Raykar PV. Heterogeneity effects on permeability-partition coefficient relationships in human stratum corneum. Pharm Res. 1988 Sep;5(9):566-573. DOI: 10.1023/A:1015989929342. PMID: 3247319.

## Equation
The framework decomposes the whole-tissue SC/water partition coefficient as a volume-weighted sum of three domains (extractable lipid, keratin protein, aqueous):

    K_SC/w = phi_lip * K_lip/w + phi_pro * PC_pro/w + phi_w

with separate power-law correlations (later recalibrated by Nitsche-Wang-Kasting 2006 from this dataset):

    K_lip/w = 0.43 * K_ow^0.81
    PC_pro/w ~ 4.2 * K_ow^0.31     (recalibrated form, used in Frasch & Barbero 2003)

The original Anderson-Raykar fit had K_lip/w ~ K_ow^~0.7 with a higher prefactor. **K-vs-Kow exponent is ~0.81 for the lipid phase and ~0.3 for the keratin/protein binding phase** -- the protein domain has much weaker lipophilicity dependence than the lipid domain. The whole-tissue K_SC/w slope vs log Kow is therefore not a single number; it is shallow at low Kow (protein-dominated) and steepens to ~0.7-0.8 at high Kow (lipid-dominated).

## Training set
The 1988-89 papers studied 21-alkyl hydrocortisone esters and other small-molecule permeants spanning approximately log Kow = -1 to +5. Subsequent Nitsche-Wang-Kasting 2006 reanalysis pooled 72 K_SC/w measurements from this and later literature.

## Reported performance
For the recalibrated (Nitsche et al. 2006) two-domain fit: rms error in log10 K_SC/w of 0.30, limited by experimental scatter.

## Validity / limitations
- Whole-tissue: hydrated, untreated human SC.
- Volume fractions for the three phases are tabulated in the original papers; phi_lip ~ 0.05, phi_pro ~ 0.45, phi_w ~ 0.5 (varies with hydration).
- Captures the curvature in K_SC/w-vs-Kow that single-power-law fits like Cleek-Bunge cannot.
- Highly polar / ionizable compounds, very large molecules (MW > 500), and metals are out of domain.

## Notes
This is the foundational decomposition that all subsequent multi-phase SC models (Mitragotri 2003, Frasch-Barbero 2003, Wang-Kasting-Nitsche 2007, Hansen 2011) build on. If you want a single equation to feed `skindiff` with a meaningful per-layer K, use the *recalibrated* coefficients from Nitsche-Wang-Kasting 2006 (separate file) rather than the 1988-89 originals -- the latter were superseded by the Klip/w = 0.43 * Kow^0.81 form on a much larger dataset.

## References
- Raykar et al. 1988 doi:10.1023/A:1015956705293
- Anderson et al. 1988 doi:10.1023/A:1015989929342
- Recalibration: Nitsche, Wang, Kasting 2006 doi:10.1002/jps.20549
