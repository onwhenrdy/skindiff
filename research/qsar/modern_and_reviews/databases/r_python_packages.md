---
name: r_python_packages
type: software (open-source libraries)
url: various — see below
license: open-source (MIT / GPL / BSD variants)
---
# Open-source R / Python packages relevant to skin permeation prediction

## Description
There is no widely-used, dedicated open-source R or Python package for skin permeation prediction prior to `skindiff` itself. Most published work is delivered as paper-supplementary scripts, ad-hoc Jupyter notebooks, or proprietary tools. The relevant ecosystem consists of general cheminformatics / ML tooling that can be combined into a skin-Kp pipeline, plus a few specific dermal-relevant packages.

## Coverage and components

### Cheminformatics descriptors
- **RDKit** (Python) — canonical open-source descriptor library (~200+ 2D, 3D descriptors plus fingerprints). https://www.rdkit.org/
- **rcdk** (R) — R bindings for the CDK Java library; ~280 2D descriptors. https://cran.r-project.org/package=rcdk
- **Mordred** (Python) — extension to RDKit with ~1800 descriptors. https://github.com/mordred-descriptor/mordred

### ML for QSAR
- **scikit-learn** (Python) — RF, SVM, GBM, ensembles. Standard reference.
- **chemprop** (Python) — message-passing graph neural networks for chemistry. https://github.com/chemprop/chemprop. Yang K et al. J Chem Inf Model 2019; 59: 3370.
- **tidymodels / caret** (R) — full ML workflow framework.

### Mechanistic skin modelling
- **`skindiff`** (R, this package) — 1D Crank-Nicolson FVM for multi-layer skin diffusion, with parameter fitting from in-vitro data.
- **deSolve / ReacTran** (R) — generic PDE solvers; can solve dermal diffusion as a custom ODE/PDE system.
- **DermalAbsorption** packages — none widely adopted as of writing; some research-code prototypes on GitHub.

### Multi-purpose QSAR / ADMET
- **OPERA** (open-source, see separate entry) — has a built-in skin Kp endpoint.
- **DeepPurpose** (Python) — deep learning for drug-target interaction; no skin endpoint by default but the framework can be repurposed.

## Access
All listed are on GitHub / CRAN / PyPI under permissive licenses (MIT, BSD, Apache, GPL). Free and open.

## Citation
Per package; see project READMEs.

## Notes
- The skin-permeation field is unusual in software terms: it has well-developed theory and reasonable QSAR datasets, but no dominant open-source reference implementation. `skindiff` aims to fill the diffusion-simulator slot. The Kp / K predictor side is still scattered across paper supplementaries.
- A practical open pipeline today: RDKit / rcdk for descriptors, scikit-learn / tidymodels for QSAR, OPERA as a sanity-check Kp predictor, `skindiff` for the mechanistic 1D simulation. All free.
