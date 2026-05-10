---
name: potts_guy_1995
predicts: Kp
units: cm/s
descriptors: MV, sum(alpha2H), sum(beta2H) (Abraham H-bond)
domain: neutral organics; small Flynn-derived set with H-bond descriptors available
status: classical
---

# Potts and Guy (1995) — Hydrogen-Bond Refinement

## Citation
Potts, R.O., Guy, R.H. (1995). A predictive algorithm for skin permeability:
the effects of molecular size and hydrogen bond activity. Pharmaceutical
Research 12(11):1628-1633. PubMed: 8592661. DOI: 10.1023/A:1016236932339.

## Equation
log Kp [cm/s] = -4.85 + 0.0256 * MV - 1.72 * sum(alpha2H) - 3.93 * sum(beta2H)

where MV is molecular volume (cm^3/mol) and sum(alpha2H), sum(beta2H) are
Abraham hydrogen-bond acidity and basicity sums.

(Coefficients as reproduced in Patel et al. 2002 Table 2 / Mitragotri 2011.)

## Training set
n ~= 37 compounds drawn from Flynn 1990 for which Abraham H-bond descriptors
were tabulated at the time. Smaller than the 1992 paper because the
descriptor coverage was the limiting factor.

## Reported performance
R^2 ~ 0.94, SD ~ 0.32 log units (small training set; better fit than the
1992 model on the same compounds). The point of the paper was mechanistic:
once H-bond descriptors are included, Kow drops out as a needed predictor.
That is, log Kow correlates with permeability *because* it correlates with
H-bond basicity, not as an independent driver.

## Validity / limitations
- Same Flynn-data caveats.
- Requires Abraham descriptors, similar to Abraham 1995/1999 -- functionally
  the same data dependency without the R2 and pi2H terms.
- Smaller dataset than other Flynn-fit models, so wider SE on coefficients.

## Notes
Mechanistically the cleaner of the two Potts-Guy papers: it argues that the
log Kow term in the 1992 form was a proxy for H-bond basicity, and by
including beta2H directly the size term cleanly captures the diffusion
contribution while the H-bond terms capture partition. Practical use is rare
because the 1992 form (only Kow + MW needed) is much easier to apply, but
this 1995 form is the conceptual predecessor of all subsequent LFER work.

## References
- Potts, R.O., Guy, R.H. (1995). Pharm. Res. 12(11):1628-1633.
  DOI: 10.1023/A:1016236932339.
- Patel, H., ten Berge, W., Cronin, M.T.D. (2002). Chemosphere 48:603-613,
  Table 2.
