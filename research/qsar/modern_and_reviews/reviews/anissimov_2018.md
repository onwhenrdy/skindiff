---
name: anissimov_2018
year: 2018
scope: review of mathematical and modelling approaches to skin permeation, with emphasis on transient and finite-dose behaviour
---
# Anissimov 2018 — review of skin permeation modelling

## Citation
Anissimov YG. Multiple review-style papers in the late 2010s — likely:
- Anissimov YG, Jepps OG, Dancik Y, Roberts MS. "Mathematical and pharmacokinetic modelling of epidermal and dermal transport processes." Advanced Drug Delivery Reviews 2013; 65: 169-190. DOI: 10.1016/j.addr.2012.11.005.
- Anissimov YG. "Mathematical models for different exposure conditions." (book chapter, ca 2018).
Verify exact citation; the user-provided "Anissimov 2018" may refer to a follow-up book chapter rather than a journal paper.

## Coverage
The 2013 ADDR paper is the more cited document and covers:
- Steady-state and transient analytical solutions to Fick's second law in single and multi-layer membranes.
- Lag-time, burst, and finite-dose corrections.
- Dispersion / variability in skin properties and how it propagates to flux uncertainty.
- Compartmental versus continuum models; when each is appropriate.
- The handling of dermal capillary clearance and the receptor-side boundary condition.
- Vehicle effects and reservoir behaviour in SC.
- Reverse problems: extracting D and K from in-vitro permeation data via non-linear fitting (relevant for `skindiff`'s `skin_fit()`).

## Key takeaways
1. Closed-form analytical solutions to single and two-layer skin permeation problems exist (Crank, Kasting) and should be used as validation benchmarks for any numerical model. `skindiff` does this in its analytical test battery.
2. Finite-dose vehicle behaviour cannot be approximated by infinite-dose Kp times duration — there's a depletion correction that grows with time.
3. Inverse fitting D and K from permeation data is ill-posed when the data is too short (steady state not reached) or too sparse; non-uniqueness is common. Penetration data (tape strip / cryocut) reduces non-uniqueness but is rarely available.
4. Stochastic / Monte Carlo simulations of skin permeation give the same mean as PDE solvers but quantify variability — useful for risk assessment.

## References
- Anissimov YG, Jepps OG, Dancik Y, Roberts MS. Adv Drug Deliv Rev 2013; 65: 169-190.
- Anissimov YG, Roberts MS. "Diffusion modelling of percutaneous absorption kinetics." J Pharm Sci 2001; 90: 504, 515; 2003; 92: 1900; 2004; 93: 470 — the foundational analytical-solution series.
