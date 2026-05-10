---
name: moss_cronin_2002
predicts: Kp
units: cm/s
descriptors: logKow, MW
domain: neutral organics; Flynn dataset with corrected steroid permeabilities
status: classical
---

# Moss and Cronin (2002) — Re-Analysis with Corrected Steroid Data

## Citation
Moss, G.P., Cronin, M.T.D. (2002). Quantitative structure-permeability
relationships for percutaneous absorption: re-analysis of steroid data.
International Journal of Pharmaceutics 238(1-2):105-109. PubMed: 11996814.
DOI: 10.1016/S0378-5173(02)00057-1.

## Equation
log Kp [cm/s] = -2.39 + 0.74 * log Kow - 0.0091 * MW

(Sign convention as published; the intercept is in cm/s after the standard
unit choice. The form is identical to Potts-Guy 1992; only the coefficients
shift because the steroid Kp values were replaced.)

## Training set
n = 116 compounds. Built on the Flynn 1990 / Wilschut 1995 compilation with
the originally-published steroid permeability data replaced by then-modern
re-measurements. Steroids had been the dominant outlier class for every
QSPR built on the Flynn set.

## Reported performance
R^2 = 0.82 (some review papers cite 0.83). Critical finding: after the
steroid replacement, residual analysis shows steroids are no longer outliers
under this model. The coefficients shift modestly versus Potts-Guy 1992
(0.74 vs 0.71 for log Kow; -0.0091 vs -0.0061 for MW; -2.39 vs the various
intercept conventions).

## Validity / limitations
- Same form as Potts-Guy, so same domain limits (neutral, aqueous donor, no
  Cleek-Bunge correction baked in).
- The "corrected" steroid Kp values are themselves debated; later
  measurements diverge again from the values Moss-Cronin used.
- Two free parameters on n=116 makes the SE on each coefficient small but
  R^2 = 0.82 leaves substantial unexplained variance.

## Notes
The most-cited "drop-in replacement" for Potts-Guy 1992 in the modern
literature: same descriptors, fitted on a cleaner dataset. If one wants a
single-equation Kp predictor and the only inputs are MW and logKow,
Moss-Cronin 2002 is arguably preferable to Potts-Guy 1992. Often quoted
side-by-side with Potts-Guy in regulatory comparisons.

## References
- Moss, G.P., Cronin, M.T.D. (2002). Int. J. Pharm. 238(1-2):105-109.
  DOI: 10.1016/S0378-5173(02)00057-1.
- Cronin et al. (1999). Investigation of the mechanism of flux across human
  skin in vitro by quantitative structure-permeability relationships.
  Eur. J. Pharm. Sci. 7:325-330. (Predecessor work.)
- Bouwman et al. (2008). Improving the applicability of (Q)SARs for
  percutaneous penetration in regulatory risk assessment. Hum. Exp. Toxicol.
  27:269-276.
