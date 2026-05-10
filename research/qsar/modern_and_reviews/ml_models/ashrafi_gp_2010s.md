---
name: ashrafi_gp
year: 2015-2020
predicts: Kp
algorithm: Gaussian Process Regression (GPR) with composite kernels
descriptors: small set of physicochemical descriptors; some studies use molecular fingerprints
training_set: extended Magnusson dataset, n ~150-300
performance: GPR test R^2 ~0.78-0.84; gives well-calibrated uncertainty estimates
---
# Ashrafi / Sun / Davey — Gaussian Process Regression for skin Kp

## Citation
Ashrafi P, Sun Y, Davey N, Adams RG, Brown MB, Moss GP. "The application of Gaussian processes in the prediction of permeability across mammalian membranes." Drug Discov Ind Pharm 2018 (verify). Earlier work in this series: Sun Y, Brown MB, Prapopoulou M, Davey N, Adams RG, Moss GP. "The application of Gaussian processes in the prediction of percutaneous absorption." J Pharm Pharmacol 2011.

## Approach
Gaussian Process Regression is well-suited to small, noisy QSAR datasets like skin Kp because it (a) does not assume a fixed functional form, (b) gives a principled uncertainty estimate at each prediction, and (c) handles non-stationary noise (different parts of descriptor space may have different prediction quality). The Sun-Ashrafi-Moss-Davey papers show GPR consistently among the top performers on Kp benchmarks. The novelty is the calibrated uncertainty: instead of just a point estimate of log Kp, the model returns a posterior distribution, which is much more useful for risk assessment than a bare number.

## Performance
Test R^2 typically 0.78-0.84 on Kp benchmarks; calibration of the predictive variance is a key advantage — coverage of nominal 95% intervals is generally close to 95%.

## Availability
GPR implementations are standard in scikit-learn, GPy, and gpytorch; the published Kp models can be reproduced from the descriptor lists in the papers.

## Limitations
GPR scales O(n^3) in training set size, so it works only for small datasets — fine for skin Kp (n < 1000) but not for combined ADMET training. Choice of kernel matters more than for tree ensembles.

## References
- Sun Y et al. J Pharm Pharmacol 2011 — original GPR for skin Kp.
- Ashrafi et al. mid-2010s follow-ups in the same series.
- Rasmussen CE, Williams CKI. Gaussian Processes for Machine Learning. MIT Press 2006 — the canonical reference.
