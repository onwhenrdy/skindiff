---
name: pawar_riviere
year: 2010-2018
predicts: Kp from realistic mixtures (vehicle effects)
algorithm: classical regression with mixture descriptors (Abraham solvation) plus PLS / SVM in later papers
descriptors: Abraham LSER solute + solvent descriptors
training_set: Riviere lab in-house porcine and human skin permeation data, hundreds of (compound, vehicle) pairs
performance: depends on paper; mixture model R^2 ~0.7-0.8 within domain
---
# Riviere lab (NCSU / Kansas State) — mixture-aware skin permeation models

## Citation
Riviere JE, Brooks JD. "Predicting skin permeability from complex chemical mixtures: dependency of quantitative structure permeation relationships on biology of skin model used." Toxicol Sci 2011; 119: 224-232. Earlier: Riviere JE, Brooks JD. Toxicol Appl Pharmacol 2005; 208: 99-110. Series continues through ~2018.

## Approach
Jim Riviere's NCSU / KSU lab is one of the few groups that systematically generated their own Kp data from porcine skin (a better human-skin surrogate than rodent skin) under matched experimental conditions, then built QSARs that explicitly include vehicle / mixture descriptors. Uses Abraham's LSER (Linear Solvation Energy Relationships) framework: solute properties (E, S, A, B, V) plus solvent properties give a five-term log-linear regression. This is the most physically grounded QSAR family for skin and the only one that handles arbitrary mixtures cleanly.

## Performance
Within the training domain (the lab's own porcine permeation database), R^2 ~0.7-0.8 with RMSE ~0.4-0.5 log units. Out-of-distribution mixtures (very different solvent compositions from training) lose predictive power — a known LSER limit.

## Availability
Data and equations in the publications; no public web tool. The Riviere lab has historically been generous with raw data on request.

## Limitations
Porcine-skin Kp is the right surrogate for some regulatory contexts but adds yet another scaling step when extrapolating to human in-vivo. Abraham descriptors require specialized software (ABSOLV / Cronnex) or hand-curation, which limits scaling to large compound libraries.

## References
- Riviere JE, Brooks JD. "Predicting skin permeability from complex chemical mixtures." Toxicol Appl Pharmacol 2005; 208: 99-110.
- Riviere JE. Chemical Mixtures and Combined Chemical and Nonchemical Stressors. CRC Press, 2018 (book chapter on dermal mixtures).
- Abraham MH. Chem Soc Rev 1993; 22: 73-83 — the original LSER paper.
