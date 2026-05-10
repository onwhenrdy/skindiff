---
name: baba_ann_2015
year: 2015
predicts: Kp
algorithm: ANN (multilayer perceptron) with feature selection
descriptors: physicochemical (MW, logP, HBD, HBA, MV, MP) plus selected topological descriptors
training_set: ~150 compounds aggregated from Flynn-derived literature
performance: R^2 around 0.80-0.85 on test set; outperforms Potts-Guy on the same compounds
---
# Baba et al. — Artificial Neural Network for Kp prediction (2015)

## Citation
Baba H, Takahara J, Mamitsuka H. "In silico predictions of human skin permeability using nonlinear quantitative structure-property relationship models." Pharmaceutical Research 2015; 32(7): 2360-2371. DOI: 10.1007/s11095-015-1629-y.

## Approach
Baba and colleagues departed from the linear Potts-Guy / Magnusson regressions by training a feed-forward neural network on a curated Flynn-derived dataset. They used a small set of physicochemical descriptors (molecular weight, logP, hydrogen-bond donors/acceptors, molar volume, melting point) augmented with topological descriptors selected by a wrapper procedure. The novelty is the explicit non-linearity: Potts-Guy assumes log Kp is a linear function of log P and MW, but skin permeation flattens at high lipophilicity (the so-called "logP plateau"). The ANN captures that without an ad hoc parabolic correction term.

## Performance
Reported test-set R^2 in the 0.80-0.85 range with leave-one-out cross-validation, compared to ~0.65-0.70 for Potts-Guy on the same data. Caveat: numbers depend on the exact data split and which compounds were excluded as outliers; the field's standard "Flynn 90 + Magnusson additions" set varies between authors.

## Availability
No public web server or repository known. The trained model is not redistributed; users would need to refit from the descriptors listed in the paper.

## Limitations
Training set still small (~150 compounds), still heavy on small organic molecules, and Kp values come from heterogeneous skin types (human cadaver, hairless mouse) and vehicle conditions. Domain of applicability is narrow.

## References
- Baba et al. 2015, Pharm Res 32: 2360-2371.
- Builds on Flynn 1990 and Magnusson et al. 2004 datasets.
