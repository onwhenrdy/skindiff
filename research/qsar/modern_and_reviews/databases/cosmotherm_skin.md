---
name: cosmotherm_skin
type: software package (commercial)
url: https://www.cosmologic.de/
license: commercial (academic licensing available at reduced cost)
---
# COSMOtherm / COSMO-RS skin partition prediction

## Description
COSMOtherm is the commercial implementation of COSMO-RS (Conductor-like Screening Model for Real Solvents), a quantum-chemical thermodynamic method for predicting partition and solvation properties from first principles. The method derives partition coefficients between any two solvents (or solvent / solid surface) from sigma profiles computed by DFT. For dermal applications, COSMO-RS-Skin or COSMOtherm-Skin gives:
- K_octanol/water (logP)
- K_SC/water (stratum corneum partition)
- K_water/lipid for various lipid classes including SC ceramides

This bypasses the need to fit empirical partition QSARs against limited data.

## Coverage
- Computes partition for any chemical with a defined molecular structure (no training-set restriction).
- Uses the "skin lipid mixture" sigma profile derived from typical SC lipid composition (ceramides, cholesterol, free fatty acids).
- Reported accuracy on K_SC/water: ~0.5-0.7 log units RMSE — comparable to or better than empirical QSARs in their training domain, and better outside it.

## Access
Commercial. Pricing for academic use is around 5-10 kEUR/year depending on configuration; industrial pricing higher. A reduced "COSMOquick" version with pre-tabulated parameters is also available.

Vendor: COSMOlogic GmbH (now part of Dassault Systèmes / BIOVIA).

## Citation
- Klamt A. "COSMO-RS: From Quantum Chemistry to Fluid Phase Thermodynamics and Drug Design." Elsevier 2005 (book).
- Klamt A, Huniar U, Spycher S, Keldenich J. "COSMOmic: a mechanistic approach to the calculation of membrane-water partition coefficients and internal distributions within membranes." J Phys Chem B 2008; 112: 12148-12157.
- Wittum R et al. (BIOVIA Cosmologic application notes for skin lipid partition).

## Notes
- One of the very few QSAR-free, mechanistic K predictors. Useful when extrapolating to compound classes outside the empirical training domain.
- Predicts only partition (K), not diffusivity (D). For D, has to be combined with an empirical D model (e.g. MW-based) or a measurement.
- Closed-source; reproducibility for a third party requires the same software and version.
- BIOVIA also markets a more recent integrated dermal absorption module that combines COSMO-RS K predictions with a 1D diffusion solver — overlaps in scope with `skindiff`.
