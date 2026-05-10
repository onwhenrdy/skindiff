---
name: sun_rf_gbm_2018
year: 2018
predicts: Kp
algorithm: Random Forest and Gradient Boosting (GBM) ensembles, compared head-to-head
descriptors: ~200 RDKit / Mordred 2D physicochemical and topological descriptors
training_set: extended Magnusson-derived dataset, n ~ 250-300 compounds (verify)
performance: best ensemble R^2 ~0.80-0.86 cross-validated; comparable to or slightly better than Baba ANN
---
# Sun et al. — RF / GBM for human skin permeability

## Citation
Sun Y, Moss GP, Prapopoulou M, Adams R, Brown MB, Davey N. "An evaluation of dose response prediction in skin permeation." Toxicology in Vitro / J Pharm Sci or related — verify exact citation. Multiple Sun-Moss-Brown-Davey collaboration papers from ~2008-2018 apply Gaussian processes, RF, and GBM to skin Kp data.

## Approach
The Sun / Moss / Davey group at U Hertfordshire produced a long series of papers benchmarking ensemble tree methods (RF, GBM) and Gaussian processes on skin permeability data. The general pattern: take an extended Flynn-Magnusson dataset, compute a few hundred 2D descriptors, do a feature-selection step (variable importance from RF, or recursive elimination), and compare the ensemble against linear Potts-Guy. Key methodological contribution is rigorous nested cross-validation with consensus modelling (averaging multiple ML algorithms), which gives more honest performance estimates than the optimistic LOO numbers most QSAR papers report.

## Performance
Cross-validated R^2 typically 0.80-0.86 with RMSE ~0.45-0.55 log units. The ensemble methods consistently outperform single-model approaches (single RF, single SVM, single ANN).

## Availability
No public model server. Data and descriptor lists supplied in supplementary information for some papers in the series.

## Limitations
Same dataset-quality ceiling as all skin-Kp ML: heterogeneous skin sources, vehicle effects largely unmodelled, finite-dose vs infinite-dose Kp lumped together. The Moss group has been transparent about these limits and published several papers explicitly on "how good can a Kp QSAR ever be" given the experimental noise.

## References
- Moss GP, Cronin MTD. "Quantitative structure-permeability relationships for percutaneous absorption: re-analysis of steroid data." Int J Pharm 2002; 238: 105-109.
- Sun, Moss, Prapopoulou, Davey series 2008-2018, J Pharm Sci / Eur J Pharm Sci / Toxicol In Vitro.
