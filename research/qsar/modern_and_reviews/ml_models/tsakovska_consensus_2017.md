---
name: tsakovska_consensus_2017
year: 2017
predicts: dermal absorption fraction; Kp via consensus of multiple models
algorithm: consensus / weighted averaging of several published QSAR models (Potts-Guy, Magnusson, Wilschut, Patel, etc.) plus mechanistic IVIVE
descriptors: standard physicochemical (MW, logP, MP, HBD/HBA)
training_set: cosmetics-relevant compound set, ~50-100 compounds with measured in-vitro Kp
performance: consensus model RMSE ~0.5 log units on test compounds; superior to any single component
---
# Tsakovska / EU CosEU consensus / mechanistic-and-QSAR hybrid

## Citation
Tsakovska I, Pajeva I, Al Sharif M, Alov P, Fioravanzo E, Kovarich S, Worth AP, Richarz AN, Yang C, Mostrag-Szlichtyng A, Cronin MTD. "Quantitative structure-skin permeability relationships." Toxicology 2017; 387: 27-42. DOI: 10.1016/j.tox.2017.06.008.

## Approach
This was a deliberate community effort under the EU SEURAT-1 / Cosmetics Europe umbrella to evaluate QSAR models for cosmetic safety assessment after the 2013 EU animal-testing ban. Rather than train a new model, the authors systematically benchmarked a dozen+ published Kp QSARs (Potts-Guy 1992, Magnusson 2004, Wilschut 1995, Patel 2002, Moss-Cronin steroid model, ten Berge 2009, etc.) on a curated cosmetics test set and assembled a consensus prediction. The mechanistic step combines QSAR-predicted Kp with a 1D diffusion model (effectively the kind of model `skindiff` implements) to get a dermal absorption fraction usable in risk assessment.

## Performance
Consensus RMSE ~0.5 log units on Kp; the consensus consistently beat any single model. The IVIVE pipeline (Kp -> diffusion model -> absorbed fraction) gave reasonable agreement with in-vivo dermal absorption studies for cosmetic ingredients.

## Availability
Models and data tables in the supplementary information of the Toxicology paper. No code release, but all component models are simple regressions and easy to reimplement.

## Limitations
Consensus only as good as its components — all built on the same noisy Flynn / Magnusson data. Domain of applicability for cosmetic ingredients (often highly lipophilic emollients, fragrances, surfactants) is shaky; many cosmetic compounds sit at the edge of the Kp QSAR training range.

## References
- Tsakovska I et al. Toxicology 2017; 387: 27-42.
- The same group also published Pajeva et al. on kinetic skin models.
- Companion EU work: Cosmetics Europe Long Range Science Strategy, ~2016-2020.
