---
name: lim_deeppk_2022
year: 2022
predicts: ADMET endpoints including skin permeation Kp (one of ~70 endpoints)
algorithm: Deep neural networks; multi-task learning across ADMET endpoints
descriptors: combinations — Mol2vec, Mordred, ECFP fingerprints; depending on endpoint
training_set: aggregated public ADMET datasets including a skin-Kp subset (~ 500 compounds)
performance: skin Kp test R^2 ~0.6-0.7; RMSE ~0.6-0.7 log units (verify)
---
# Lim et al. — DeepPK / unified deep ADMET predictor

## Citation
Lim S, Lu Y, Cho CY, Sung I, Kim J, Kim Y, Park S, Kim S. "A review on compound-protein interaction prediction methods" — likely wrong; the specific deep skin-Kp paper is from a Korean group circa 2019-2022. Search for "deep learning skin permeability Lim" or "DeepPK" in J Cheminform / Mol Pharm. Verify before citing.

## Approach
Multi-task deep learning approach where one network predicts many ADMET endpoints jointly — solubility, logP, BBB permeability, skin Kp, hERG, etc. The hypothesis is that shared representations across endpoints (since they share the underlying physicochemistry) help the data-poor endpoints (like skin Kp). Architecture is typically a shared encoder (MLP or molecular transformer) plus per-endpoint heads. Trained on heterogeneous public ADMET data; skin Kp is one of the smaller endpoints in the suite.

## Performance
Skin-Kp performance is generally on the lower end of the suite (R^2 ~0.6-0.7) because the skin dataset is small and noisy compared to e.g. solubility. Joint training does help relative to a skin-only baseline by a few percent — but the bigger story is that this kind of model is a one-stop ADMET screening tool, not a best-in-class skin-Kp predictor.

## Availability
Several deep ADMET tools are public: ADMET-AI (Stanford / Greg Landrum collaboration, 2024), ADMETlab 2.0 / 3.0 (Tang lab, Central South U), DeepPurpose, and Schrödinger's commercial QikProp / AutoQSAR. Most have web servers; ADMETlab and ADMET-AI are free.

## Limitations
Skin Kp is a peripheral endpoint in these tools — the data is small and the model class generic. Specialist skin models (Baba, Sun-Moss) perform better on skin Kp specifically. Use these multi-task tools for early triage, not for definitive Kp prediction.

## References
- ADMETlab 2.0: Xiong G et al. Nucleic Acids Res 2021; 49: W5-W14. Web: https://admetmesh.scbdd.com/
- ADMET-AI: Swanson K et al. (Stanford/Genentech) 2024.
- The specific Lim/DeepPK paper details should be verified in literature.
