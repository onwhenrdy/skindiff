---
name: potts_guy_1992_implied_K
predicts: K (implicit, derived from K_p decomposition)
phase_pair: SC/water (whole tissue, implicit)
units: dimensionless
descriptors: logKow, MW
domain: ~93 compounds in the Flynn 1990 dermal absorption database; log Kow -3 to +6, MW 18-750 Da
status: classical (the foundational K_p QSAR; K is implicit, not directly fit)
---

# Potts & Guy (1992) -- implicit K_SC/water from skin permeability

## Citation
Potts RO, Guy RH. Predicting skin permeability. Pharm Res. 1992 May;9(5):663-669. DOI: 10.1023/A:1015810312465. PMID: 1608900.

## Equation
Potts and Guy fit the in-vitro skin permeability coefficient K_p as

    log K_p (cm/s) = -2.74 + 0.71 * log K_ow - 0.0061 * MW

If permeability is decomposed as K_p = K_SC/w * D_SC / h_SC and one assumes the diffusion coefficient D_SC depends only on MW (free-volume picture: log D ~ -alpha * MW), then the *implicit* SC/water partition coefficient slope vs log Kow is the slope on log K_p, namely

    K_SC/w  ~  K_ow^0.71

so the **implicit K-vs-Kow exponent in the Potts-Guy K_p model is ~0.71**. Note this is *implicit* -- Potts and Guy did not separately measure or fit K_SC/w; the 0.71 slope is what a reader infers when separating the K and D contributions. Subsequent work (Cleek-Bunge, Anderson-Raykar, Wang-Kasting) found 0.74-0.81 when fitting K_SC/w directly.

## Training set
Approximately 93 compounds from the Flynn 1990 review of in-vitro human and mammalian skin permeability data. log Kow range -3 to +6, MW range 18 to >750 Da. Mostly small organic permeants from aqueous donor solutions.

## Reported performance
On log K_p: R^2 ~ 0.67, RMSE ~ 0.7 in log10 K_p across the Flynn 1990 set. K_p is the directly fit quantity; performance on K_SC/w is not separately reported (because K_SC/w is implicit).

## Validity / limitations
- K_p (not K) is the fit quantity; readers who want a K from this paper are doing arithmetic, not citing a measurement.
- log Kow [-3, +6], MW < ~ 750 Da, neutral organics, aqueous donor.
- Single-power-law in Kow (no curvature for very lipophilic / hydrophilic regimes), no decomposition of lipid vs corneocyte phases.

## Notes
The Potts-Guy paper is the founding QSAR for transdermal absorption in general. For K *specifically*, it is a useful cross-check: the implicit K_SC/w slope of 0.71 sits near the consensus middle of the literature (0.67 Hansen, 0.69 Frasch-Barbero, 0.74 Cleek-Bunge, 0.81 Wang-Kasting). For `skindiff`, do not use Potts-Guy directly to set a K -- it conflates K with D. Use Cleek-Bunge if you want a single-power-law K_SC/w from log Kow alone, or Nitsche-Wang-Kasting / Hansen for a two-domain K_SC/w.

## References
- Potts & Guy 1992 doi:10.1023/A:1015810312465
- Flynn 1990 (parent dataset): GL Flynn, in *Principles of Route-to-Route Extrapolation for Risk Assessment*, Elsevier (1990) pp 93-127.
