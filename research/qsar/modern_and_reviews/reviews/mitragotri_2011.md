---
name: mitragotri_2011
year: 2011
scope: comprehensive review of mathematical models of skin permeability — from steady-state Kp regressions to compartmental and finite-element diffusion models
---
# Mitragotri et al. 2011 — Mathematical models of skin permeability

## Citation
Mitragotri S, Anissimov YG, Bunge AL, Frasch HF, Guy RH, Hadgraft J, Kasting GB, Lane ME, Roberts MS. "Mathematical models of skin permeability: an overview." International Journal of Pharmaceutics 2011; 418: 115-129. DOI: 10.1016/j.ijpharm.2011.02.023.

## Coverage
Co-authored by basically every senior figure in the skin-modelling field as a community consensus document. Covers:
- Steady-state and transient permeation theory (Fick's laws, lag-time, K_p / D / K relationships).
- Empirical Kp QSARs (Potts-Guy, Magnusson, etc.) — context and limitations.
- Compartmental pharmacokinetic models (vehicle / SC / viable epidermis / dermis / receptor).
- Finite-difference and finite-element 1D diffusion models — the kind `skindiff` implements.
- Multi-layer, multi-pathway models (lipid lamellae, corneocyte intracellular, shunt routes).
- Stratum corneum heterogeneity (brick-and-mortar geometry, lipid pathway tortuosity).
- In vitro / in vivo comparisons and the IVIVE problem.
- Vehicle effects and finite-dose behaviour.

## Key takeaways
1. The field is methodologically mature: there is no shortage of models. The shortage is in well-controlled data to validate them against.
2. Compartmental pharmacokinetic models and 1D PDE diffusion models give equivalent steady-state predictions but differ in lag-time and finite-dose dynamics; for risk assessment with finite-dose exposures, the 1D PDE approach is the right tool.
3. Microscopic brick-and-mortar models add interpretability but rarely change the macroscopic prediction enough to justify the added complexity unless the question is specifically about the SC microstructure.
4. The lipid pathway dominates Kp for lipophilic compounds; for very polar compounds the polar pathway and shunt route (hair follicles) become non-negligible.
5. QSAR Kp + mechanistic diffusion model is the recommended pipeline for risk assessment, with Kp providing the boundary condition.

## References
The paper itself cites ~250 references. Key foundational works:
- Crank J. The Mathematics of Diffusion, 2nd ed. Oxford 1975 — the standard PDE reference.
- Kasting GB. J Pharm Sci 2001; 90: 202 — finite-dose membrane permeation analytical solutions (validated against in `skindiff`'s test suite).
- Anissimov YG, Roberts MS series on skin permeation theory, J Pharm Sci 2001-2010.
- Frasch HF, Barbero AM. J Pharm Sci 2008 — random-walk SC model.
