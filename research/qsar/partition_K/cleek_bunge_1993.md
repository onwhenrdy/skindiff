---
name: cleek_bunge_1993
predicts: K
phase_pair: SC/water (whole tissue)
units: dimensionless (mass-based concentration ratio at equilibrium; original equations sometimes presented in mL/g for SC mass)
descriptors: logKow
domain: small organic permeants used in the Flynn 1990 dermal absorption database; broadly log Kow approximately -3 to +6, MW < 500 Da
status: classical
---

# Cleek & Bunge (1993) -- single-power-law SC/water partition

## Citation
- Cleek RL, Bunge AL. A new method for estimating dermal absorption from chemical exposure. 1. General approach. Pharm Res. 1993 Apr;10(4):497-506. DOI: 10.1023/A:1018981515480.
- Bunge AL, Cleek RL. A new method for estimating dermal absorption from chemical exposure. 2. Effect of molecular weight and octanol-water partitioning. Pharm Res. 1995 Jan;12(1):88-95. DOI: 10.1023/A:1016242821610. PMID: 7724493.

## Equation
The Cleek-Bunge default (paired with the Potts-Guy 1992 permeability correlation) for the SC/water partition coefficient is the simple single power law:

    log K_SC/w = b * log K_ow + a
    with a = 0, b = 0.74    (-->  K_SC/w = K_ow^0.74)

so the **K-vs-Kow exponent is 0.74**. (Compare: implicit Potts-Guy 1992 exponent ~0.71, Nitsche/Wang/Kasting 2006 lipid-only exponent 0.81, Hansen/Frasch-Barbero lipid exponent 0.69-0.70.)

This is implemented as a one-line update to a Potts-Guy K_p estimate via a B-factor correction that re-introduces a viable-epidermis resistance for very lipophilic compounds; the K_SC/w piece itself is the simple log-log line above.

## Training set
Calibrated against the Flynn 1990 in-vitro permeability database (~90+ compounds collated by Flynn for the original Potts-Guy regression).

## Reported performance
RMS error in log10 K_SC/w of order 0.4-0.6, characteristic of single-power-law fits to a heterogeneous tissue. Predictions degrade systematically for very lipophilic (log Kow > 4) and very hydrophilic (log Kow < -1) compounds because the corneocyte / aqueous contributions are absorbed into a single slope.

## Validity / limitations
- MW range: tested up to ~750 Da, recommended for MW < 500 Da.
- log Kow range: empirically -3 to +6 in the training set; predictions for log Kow outside ~[-1, 4] are unreliable because the underlying K_SC/w is two-domain in nature and a single line cannot capture the curvature.
- Whole-tissue (hydrated SC); does not separate lipid and corneocyte phases.
- Charged species, large biomolecules, metals: out of domain.

## Notes
Cleek-Bunge is the workhorse "default" K_SC/w correlation in regulatory and risk-assessment dermal-absorption tools (US EPA dermal-exposure guidance, IH SkinPerm-style spreadsheets) precisely because it is one number. For mechanistic transport models like `skindiff` that allow a per-layer K, the two-domain Anderson-Raykar / Nitsche-Wang-Kasting decomposition is preferable -- but Cleek-Bunge is the easy one-line K to start with when only log Kow is known.

## References
- Cleek & Bunge 1993 doi:10.1023/A:1018981515480
- Bunge & Cleek 1995 doi:10.1023/A:1016242821610
- Potts & Guy 1992 (parent K_p correlation) doi:10.1023/A:1015810312465
