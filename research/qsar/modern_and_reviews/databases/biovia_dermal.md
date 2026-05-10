---
name: biovia_dermal
type: software package (commercial)
url: https://www.3ds.com/products-services/biovia/products/scientific-informatics/cosmologic/
license: commercial
---
# BIOVIA Cosmologic dermal absorption modules

## Description
After Dassault Systèmes acquired COSMOlogic GmbH (2020), the COSMO-RS / COSMOtherm portfolio was integrated into BIOVIA's product line. The dermal-absorption-relevant components include:
- **COSMO-RS-Skin** — sigma-profile-based prediction of skin lipid / corneocyte partition coefficients.
- **COSMOmic** — micelle / membrane partition prediction; can model SC lipid lamellae.
- **COSMOplex** — newer (2018+) extension for self-organizing systems.
- A higher-level dermal absorption module that combines COSMO-RS K predictions with a 1D diffusion solver to predict dermal absorption fractions over time. Functionally similar to what `skindiff` does but with COSMO-RS K predictions feeding the diffusion model.

Also from BIOVIA: ADMET Predictor / Pipeline Pilot ADMET nodes (different lineage, ex-Accelrys), with their own skin permeation endpoints.

## Coverage
- Mechanistic K prediction for arbitrary chemicals — no QSAR training-set restriction.
- Includes lipid / water / corneocyte partition predictions tuned to skin composition.
- Diffusion module covers single-vehicle finite-dose, multi-vehicle, and occluded scenarios.

## Access
- Commercial; pricing on application. Academic licensing through Dassault Systèmes' academic program.
- Closed source. Integrates with Materials Studio / Pipeline Pilot ecosystems.
- Web hosted "BIOVIA Cosmocloud" option for pay-per-use.

## Citation
- Klamt A, Eckert F, Arlt W. "COSMO-RS: An Alternative to Simulation for Calculating Thermodynamic Properties of Liquid Mixtures." Annu Rev Chem Biomol Eng 2010; 1: 101.
- BIOVIA application notes on skin lipid partition (downloadable from 3ds.com).

## Notes
- The COSMO-RS-Skin K predictor is one of the few mechanistic alternatives to empirical Kp/K QSARs and has the unique advantage of being applicable to compounds far outside any empirical training set.
- The integrated dermal-absorption module overlaps in scope with `skindiff` (Cosmologic does the 1D diffusion solve too) but with closed-source K prediction.
- For users without commercial licensing, `skindiff` + an open Kp/K source (OPERA, published QSAR equations) covers similar territory.
