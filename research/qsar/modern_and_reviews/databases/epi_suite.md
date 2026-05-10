---
name: epi_suite
type: software package (free, closed source)
url: https://www.epa.gov/tsca-screening-tools/epi-suitetm-estimation-program-interface
license: free, U.S. EPA-provided; closed source
---
# EPI Suite — KOWWIN, DERMWIN, and friends

## Description
EPI Suite (Estimation Programs Interface) is the long-running U.S. EPA chemical-property estimation suite, originally developed by Syracuse Research Corporation (SRC) and maintained by EPA. It contains ~10 standalone modules for different physicochemical and environmental fate endpoints. The two relevant for dermal work are:

- **KOWWIN**: log octanol-water partition coefficient (logP) prediction — fragment-additivity method (Hansch-Leo / Meylan). Among the most cited logP predictors and a near-universal default.
- **DERMWIN**: dermal absorption-rate prediction. Uses the Robinson-modified-Potts-Guy equation (a logP and MW based skin Kp QSAR) and an exposure scenario calculator. Outputs dermal Kp and a 1-hour absorption estimate per cm^2.

## Coverage
- KOWWIN: trained on ~10000+ compounds; broad applicability.
- DERMWIN: built around the Potts-Guy / Robinson Kp equation; doesn't have its own training data per se, just applies the published equation to user-supplied SMILES.
- Other modules (BIOWIN, AOPWIN, BCFBAF) cover environmental fate and are not dermal-specific.

## Access
- Free download from EPA: https://www.epa.gov/tsca-screening-tools/epi-suitetm-estimation-program-interface
- Windows desktop application; runs under Wine on Linux/macOS.
- No API or library — interactive GUI plus batch-input mode via text files.
- Closed source; cannot be embedded in other tools.

## Citation
US EPA. EPI Suite v4.11. Estimation Programs Interface Suite for Microsoft Windows. United States Environmental Protection Agency, Washington, DC, USA. 2012.

For DERMWIN specifically: Robinson PJ. "Methods for the Estimation of Dermal Absorption: Recommendations for ATSDR." 1996. (DERMWIN's underlying method.)

## Notes
- DERMWIN's Kp is just Potts-Guy with a few modifications - functionally equivalent to running the equation by hand.
- KOWWIN logP predictions are standard inputs for downstream Kp QSARs and are reasonably accurate for most drug-like molecules (RMSE ~0.3-0.4 log units).
- Maintained and supported by EPA; updates are infrequent but reliable.
- For modern open-source workflow, OPERA (above) supersedes EPI Suite for most tasks while remaining free and openly licensed.
