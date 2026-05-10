---
name: heeds_dermal
type: software package (commercial)
url: https://www.3ds.com/products-services/simulia/products/isight-simulia-execution-engine/heeds/
license: commercial
---
# HEEDS (Dassault Systèmes) — design exploration; not a dermal-specific tool

## Description
HEEDS is a Dassault Systèmes design-exploration / multidisciplinary optimization platform — it wraps around external simulation tools and runs parameter sweeps, design-of-experiments, and optimization. The user-provided "HEEDS dermal absorption" is most likely a reference to HEEDS being used to drive a third-party dermal absorption model (e.g. an Abaqus / SIMULIA dermal FEA simulation, or a coupled COSMOlogic + diffusion model) rather than HEEDS having native dermal capability.

If the user has seen HEEDS marketed as a dermal tool, it would be in the context of a packaged SIMULIA / BIOVIA workflow — HEEDS as the optimization engine, BIOVIA / COSMOlogic as the property predictor, and a simulation tool (SIMULIA, Abaqus, OpenFOAM) as the diffusion solver. Verify the specific product and lineage with Dassault Systèmes if this is being seriously considered.

## Coverage
- Generic — depends entirely on what underlying simulation tool HEEDS is driving.
- Not a Kp / K predictor itself.

## Access
- Commercial; SIMULIA pricing.
- Strong scripting / API access (Python, batch).

## Citation
Dassault Systèmes SIMULIA Corp. *HEEDS Documentation*. https://www.3ds.com/

## Notes
- This entry is included for completeness given the user's source list, but HEEDS is a parameter-sweep / optimization engine, not a dermal model. If you actually need optimization over a dermal-absorption simulation, `skindiff`'s `skin_fit()` plus standard optimization / sensitivity packages in R (DEoptim, sensitivity, mco) cover the same ground for free.
- The Dassault dermal story is properly the BIOVIA / COSMOlogic stack (separate entry); HEEDS sits on top of it for design exploration.
