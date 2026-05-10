---
name: hansen_2011_keratin
predicts: K
phase_pair: corneocyte/water (specifically keratin/water binding, PC_pro/w in two-domain framework)
units: dimensionless (mass-based -- K_Nernst defined as bound mass per dry-keratin mass divided by free aqueous concentration; can be converted to volume-based using corneocyte protein density)
descriptors: logD (logKow corrected for ionization at pH 7.4)
domain: 64 compounds spanning log D -4.67 to +5.70 and MW 32-1373 Da
status: modern (benchmark for keratin-binding QSAR)
---

# Hansen, Selzer, Schaefer & Kasting (2011) -- extended keratin-binding database

## Citation
Hansen S, Selzer D, Schaefer UF, Kasting GB. An extended database of keratin binding. J Pharm Sci. 2011 May;100(5):1712-1726. DOI: 10.1002/jps.22396. PMID: 21246561.

## Equation
Linear free-energy relationship for keratin-water binding (covering hoof-and-horn powder, callus, delipidized SC):

    log K_Nernst = 0.65 (+/- 0.07) + 0.32 (+/- 0.02) * log D

where D is the pH-corrected octanol-water partition coefficient (logD at pH 7.4, equivalent to logKow for fully neutral compounds).

So **K-vs-Kow exponent for keratin/protein binding is 0.32** -- very close to the Nitsche-Wang-Kasting protein-phase exponent of 0.31 and the Frasch-Barbero protein exponent of 0.31. The keratin/protein domain is much *less* lipophilicity-sensitive than the lipid domain (slopes 0.67-0.81 across the literature for lipid).

## Training set
64 chemically diverse compounds in the unified, error-corrected dataset (75 in the original "extended database" before curation). MW 32-1373 Da. log D -4.67 to +5.70. Sources include hoof-and-horn powder binding studies, callus uptake studies, and delipidized SC partition measurements.

## Reported performance
R^2 = 0.78, standard error SE = 0.35 in log K_Nernst.

## Validity / limitations
- K_Nernst is mass-based (mass bound per mass dry keratin) and must be converted with phi_pro / rho_pro before substitution into a volume-based whole-tissue K_SC/w sum.
- Some compounds bind keratin via specific mechanisms not captured by lipophilicity alone (azoles, metal-binders, etc.) and these are outliers in the fit.
- The small slope (0.32) means: for moderately lipophilic compounds the protein domain is typically a smaller K-contribution than the lipid domain; for hydrophilic compounds (log D < 0) it is dominant alongside the aqueous phase.

## Notes
This is the canonical recalibration of the protein/keratin K-side correlation. Together with the Hansen 2013 lipid correlation (K_lip/w = K_ow^0.67) it gives a complete two-domain K_SC/w predictor with internally consistent calibration. For `skindiff`, this is the value to use for a corneocyte-layer K when the corneocyte is modeled as a separate compartment.

## References
- Hansen, Selzer, Schaefer, Kasting 2011 doi:10.1002/jps.22396
- Lipid companion: Hansen, Lehr, Schaefer 2013 doi:10.1016/j.addr.2012.04.011
