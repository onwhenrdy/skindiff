---
name: lien_gao_1995
predicts: Kp
units: cm/s
descriptors: logKow, (logKow)^2, MW (or log MW), Hb (hydrogen-bond count)
domain: human and hairless-mouse skin, 53 miscellaneous compounds, aqueous vehicle
status: classical
---

# Lien and Gao (1995) — Bilinear Kow + Hydrogen-Bond Term

## Citation
Lien, E.J., Gao, H. (1995). QSAR analysis of skin permeability of various
drugs in man as compared to in vivo and in vitro studies in rodents.
Pharmaceutical Research 12(4):583-587. PubMed: 7596996.
DOI: 10.1023/A:1016266316100.

## Equation
log Kp [cm/s] = 0.84 * log Kow - 0.07 * (log Kow)^2 - 0.27 * Hb - 1.84 * log MW + 0.0833

(Coefficients as quoted in the PMC review by Wang et al. 2018 / Mitragotri
2011 reproductions; the original paper gives the form
log P_perm = a*logKow + b*(logKow)^2 - c*Hb - d*logMW + e with `a, b, c, d`
all positive.)

## Training set
n = 53 compounds spanning a broad chemical-class miscellany; in-vitro
permeability through human skin (and a parallel hairless-mouse fit). The
authors built the same regression on rodent data to allow cross-species
extrapolation.

## Reported performance
Original paper R values around 0.92-0.95 for the human-skin fit (the model
has 4 free parameters on n=53, so R^2 is generous). Adding (logKow)^2 and an
explicit hydrogen-bond count was the major change vs Potts-Guy: it captures
the often-observed flattening / inversion of Kp at very high logKow.

## Validity / limitations
- Hb is a simple count of donor + acceptor groups; substantially less rich
  than Abraham's A and B descriptors but easier to compute by hand.
- (logKow)^2 term means the model has a *maximum* in Kow space (around
  log Kow ~ 6), beyond which it predicts decreasing Kp -- this matches the
  Cleek-Bunge intuition that very lipophilic compounds are limited by
  aqueous-epidermis resistance, but is empirical.
- 53 compounds is small; the cross-species fit is weaker.
- Steroid outliers persist.

## Notes
Lien-Gao is one of the first papers to add a hydrogen-bond term to the basic
Potts-Guy form. It was followed by Patel-ten Berge-Cronin 2002 (which also
includes hydrogen bonding), Abraham 1995/1999/2004 (full LFER), and
Moss-Cronin 2002 (steroid re-fit). Useful as an interpretable "next step
beyond Potts-Guy" but largely superseded by Abraham's LFER for accuracy and
by ten Berge SKINPERM for engineering use.

## References
- Lien, E.J., Gao, H. (1995). QSAR analysis of skin permeability of various
  drugs in man as compared to in vivo and in vitro studies in rodents.
  Pharm. Res. 12(4):583-587. DOI: 10.1023/A:1016266316100.
- Reproductions of the equation: Mitragotri et al. (2011). Int. J. Pharm.
  418:115-129. DOI: 10.1016/j.ijpharm.2011.02.023.
- Lien, E.J., Tong, G.L. (1973). Physicochemical properties and percutaneous
  absorption of drugs. J. Soc. Cosmet. Chem. 24:371-384. (Predecessor work.)
