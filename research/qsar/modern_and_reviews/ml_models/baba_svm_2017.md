---
name: baba_svm_2017
year: 2017
predicts: Kp
algorithm: SVM (support vector regression) with RBF kernel; ensembled
descriptors: 2D MOE descriptors with feature selection
training_set: aggregated Flynn + Magnusson + literature additions, on the order of 200-250 compounds
performance: test-set R^2 ~ 0.83; RMSE ~0.5 log units
---
# Baba et al. — SVM ensembles for Kp (2017)

## Citation
Baba H, Takahara J, Yamashita F, Hashida M. "Modeling and prediction of solvent effect on human skin permeability using support vector regression and random forest." Pharmaceutical Research 2017; 34(9): 1854-1865. DOI: 10.1007/s11095-017-2197-0 (verify DOI; details from author group's series of skin-Kp ML papers).

## Approach
A follow-up to the 2015 ANN work that switched the learner to support-vector regression and bagged ensembles, and crucially attempted to model the *vehicle effect* — i.e. how Kp changes when the same compound is dosed from water vs. ethanol vs. propylene glycol. Most QSAR models assume "infinite-dose aqueous donor", which is wildly unrepresentative of real cosmetic / topical drug exposures. Baba added solvent descriptors (Hansen solubility parameters, dielectric constant) alongside the solute descriptors so the model can take a (compound, solvent) pair as input.

## Performance
Test-set R^2 on the order of 0.83, RMSE near 0.5 log units. Performance was robust to leave-one-vehicle-out cross-validation, suggesting some genuine vehicle-effect generalization. Numbers should be checked against the paper.

## Availability
No public model release known. The descriptor lists and training data structure are documented in the paper.

## Limitations
Vehicle coverage still skewed toward the common lab solvents (water, ethanol, PG, mineral oil). Skin source (human vs. animal) heterogeneity not corrected for. Like all skin-Kp models, struggles with very lipophilic (logP > 5) or very polar (logP < -1) compounds — both regions are sparsely populated in training.

## References
- Baba et al. 2017, Pharm Res — verify exact issue/DOI.
- Earlier ANN paper: Baba et al. 2015.
- Vehicle-effect rationale: Twist & Zatz 1986; Cross & Roberts 1995.
