---
name: logKow_predictors_overview
predicts: K (octanol/water) -- input to all the SC K QSARs above
phase_pair: octanol/water
units: dimensionless (log10 of mole-fraction-based or volume-based K_ow; conventionally reported as logP / log Kow)
descriptors: full molecular structure
status: standard reference inputs
---

# logKow predictors (input to the skin K models)

The skin-K QSARs in this folder all take log Kow as their primary descriptor. When experimental log Kow is unavailable, predicted values are routinely substituted -- but the predictors disagree by 0.3-1.0 log units on hard cases, and that disagreement propagates into K_SC/water through the slope (~0.7), so the choice of logKow predictor matters.

## Common predictors

### KOWWIN (US EPA EPI Suite)
- Atom/fragment contribution method, fragment-based QSAR.
- Trained on ~13 000 compounds from the SRC PHYSPROP database.
- Reported performance: ~67% of predictions within +/- 1 log unit of experiment; ionisable acids (slope ~1.24) and bases (slope ~0.66) are systematically biased.
- Distributed as part of EPA EPI Suite (free, Windows installer).
- Reference: https://www.epa.gov/tsca-screening-tools/epi-suite-tm-estimation-program-interface

### ALOGPS (Tetko / VCCLAB)
- Neural-network ensemble trained on ~12 000 PHYSPROP compounds plus newer literature.
- Reported RMSE ~0.4-0.5 log units on independent test sets.
- Available online at vcclab.org and as a downloadable client.
- Reference: Tetko et al. 2005 J Chem Inf Model 45:1382. doi:10.1021/ci049854h

### ClogP (Biobyte / Hansch & Leo)
- Fragment-based with proximity corrections; the historical gold standard fragment method.
- Trained on the Pomona/Biobyte MedChem database (~10 000+ compounds).
- Reported RMSE ~0.4 log units on benchmark sets; tends to be the most accurate for typical drug-like molecules.
- Commercial; available as standalone software and within ChemDraw.
- Reference: Leo AJ. Chem Rev. 1993;93:1281. doi:10.1021/cr00020a002

### miLogP (Molinspiration)
- Group contribution method calibrated on a curated drug-like set of ~12 000 molecules.
- Free web service at molinspiration.com.
- Reported RMSE ~0.5-0.7 log units; intended as a fast first estimate, not a precision tool.
- Reference: https://www.molinspiration.com/services/logp.html

### XLogP3 (Cheng et al.)
- Atom-additive method with knowledge-based corrections.
- Trained on ~8 000 reliable measurements; benchmark RMSE ~0.7 log units on FDA-approved drugs (independent set).
- Free; integrated into PubChem so every public PubChem CID has an XLogP3 value.
- Reference: Cheng T et al. 2007 J Chem Inf Model 47:2140. doi:10.1021/ci700257y

## Notes for `skindiff`

- For `skindiff`, the practical pipeline is: experimental log Kow -> use it; else default to whichever predictor has been validated on your compound class (drug-likes: ClogP or XLogP3; environmental chemicals: KOWWIN; surfactants / amphiphiles: be very cautious -- single-Kow predictors break down).
- All these predictors return a single number -- they do not give an uncertainty. For an honest skindiff K input, propagate at minimum +/- 0.5 log units of Kow uncertainty into the K used for forward simulation.
- Ionizable compounds: use the *neutral* logKow (not logD at the relevant pH) and apply the partitioning correction at the K_SC/w step, *unless* the K-side QSAR explicitly takes log D as input (e.g., Hansen 2011 keratin model uses log D at pH 7.4).
- These are octanol/water predictors, not skin-K predictors. Always pair them with a skin-K QSAR (Cleek-Bunge / Hansen / Wang-Kasting / Anderson-Raykar) to get the K used by `skindiff`.

## References
- Mannhold R, Poda GI, Ostermann C, Tetko IV. Calculation of molecular lipophilicity: state-of-the-art and comparison of log P methods on more than 96,000 compounds. J Pharm Sci. 2009 Mar;98(3):861-93. doi:10.1002/jps.21494.
- Plante J, Werner S. JPlogP: an improved logP predictor trained using predicted data. J Cheminform. 2018 Dec 6;10(1):61. doi:10.1186/s13321-018-0316-5.
