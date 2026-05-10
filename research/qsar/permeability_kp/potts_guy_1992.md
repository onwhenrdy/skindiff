---
name: potts_guy_1992
predicts: Kp
units: cm/s
descriptors: logKow, MW
domain: stratum corneum, neutral organics, MW 18-750, logKow -3 to +6, aqueous vehicle
status: classical
---

# Potts and Guy (1992) — Predicting Skin Permeability

## Citation
Potts, R.O. and Guy, R.H. (1992). Predicting skin permeability. Pharmaceutical Research 9(5):663-669. DOI: 10.1023/A:1015810312465. PubMed: 1608900.

## Equation
log Kp [cm/s] = -6.3 + 0.71 * log Kow - 0.0061 * MW

(Sometimes also written as log Kp [cm/s] = 0.71 * log Kow - 0.0061 * MW - 2.7 in
the c-formulation; the -6.3 form is the standard log10 cm/s intercept after
combining log(D0/h) and the activity scaling. The -2.7 vs -6.3 split is
discussed in Mitragotri 2011.)

## Training set
n = 93 compounds (paper says "more than 90") drawn from Flynn (1990)'s
compilation of in-vitro permeation through human skin from aqueous solution.
MW range: 18 to >750 Da. log Kow range: -3 to +6. The Flynn 1990 dataset is
the foundational benchmark and has been re-analyzed dozens of times.

## Reported performance
R^2 = 0.67. The remaining ~33% variance has driven essentially every
follow-up paper in this list. Authors did not report cross-validated Q^2 or
prediction-set RMSE in the original.

## Validity / limitations
- Aqueous vehicle assumed. Not directly applicable to organic solvents,
  semisolids, or vehicles that perturb the SC.
- Treats SC as a single homogeneous lipid membrane. No aqueous-boundary-layer
  correction; for very lipophilic solutes (log Kow > ~4) the model
  over-predicts Kp because the viable epidermis becomes rate-limiting -- this
  is the gap Cleek-Bunge 1993 fills.
- Neutral, non-ionizable molecules. No pH/ionization handling.
- Steroid sub-class is a known outlier population (re-analyzed by Moss-Cronin 2002).
- MW is a coarse proxy for molecular volume; for halogenated compounds with
  high density, Vecchia-Bunge 2002 recommends adjusting MW by liquid density.

## Notes
The reference QSAR for skin permeability. Adopted by US EPA (1992 Dermal
Exposure Assessment), NIOSH (Skin Permeation Calculator), and almost every
risk-assessment tool. The model says SC intercellular lipid properties alone
suffice to capture Kp's MW and Kow dependence -- a result with major
mechanistic implications. R^2 = 0.67 is the bar every later model is
benchmarked against.

## References
- Potts, R.O., Guy, R.H. (1992). Predicting skin permeability. Pharm. Res.
  9(5):663-669. DOI: 10.1023/A:1015810312465.
- Flynn, G.L. (1990). Physicochemical determinants of skin absorption. In:
  Gerrity, T.R., Henry, C.J. (Eds.), Principles of Route-to-Route
  Extrapolation for Risk Assessment. Elsevier, New York, pp. 93-127.
- Mitragotri et al. (2011). Mathematical models of skin permeability: an
  overview. Int. J. Pharm. 418:115-129. DOI: 10.1016/j.ijpharm.2011.02.023.
- VEGA QMRF: https://www.vegahub.eu/vegahub-dwn/qmrf/QMRF_SKIN_KP_POTTS.pdf
