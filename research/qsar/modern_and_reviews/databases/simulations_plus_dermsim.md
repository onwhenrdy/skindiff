---
name: simulations_plus_dermsim
type: software package (commercial)
url: https://www.simulations-plus.com/
license: commercial
---
# Simulations Plus DermSim / GastroPlus dermal module

## Description
Simulations Plus is a U.S. pharmaceutical-modelling software vendor whose flagship product, GastroPlus, is the industry-standard PBPK simulator for oral absorption. They also offer:
- **DermSim** (or transdermal module of GastroPlus, depending on configuration): a mechanistic dermal absorption simulator. Models the vehicle, stratum corneum, viable epidermis, dermis, and capillary uptake in a 1D compartmental scheme similar to what `skindiff` does, with built-in QSAR predictions for D, K, and Kp.
- **ADMET Predictor** (separate product, formerly ex-Pharsight): high-throughput ADMET QSAR engine with skin Kp, skin sensitization, and other dermal endpoints.

## Coverage
- DermSim: any compound the user can describe in terms of physicochemical properties; built-in QSARs for D and K from the structure.
- ADMET Predictor: ~150 endpoints including skin Kp; trained on extensive proprietary + public data.
- Used widely in pharmaceutical and FMCG industry; FDA familiar with GastroPlus dermal output for transdermal patch submissions.

## Access
- Commercial licensing through Simulations Plus.
- Closed source. Strong professional support, training, and validation packages.
- API and scripting available within their software environment; less convenient than open libraries for embedding in custom pipelines.

## Citation
- Simulations Plus, Inc. *GastroPlus / DermSim User Manual*. Lancaster, CA, USA. https://www.simulations-plus.com/
- Lukacova V et al. "Approximation of the partition coefficient of solutes between aqueous and lipid phases of the stratum corneum using ADMET Predictor." (Conference proceedings; verify exact citation.)

## Notes
- Closest commercial peer to `skindiff`'s mechanistic simulator role, with the addition of bundled QSAR for D and K and bundled PBPK for systemic translation.
- Validation packages for FDA submissions are a major value-add for regulated work.
- For research / academic / open-pipeline contexts, `skindiff` provides the diffusion-simulator equivalent under an open license, leaving the user free to choose any QSAR for D and K (open or otherwise).
