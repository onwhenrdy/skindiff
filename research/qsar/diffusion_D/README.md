# Diffusion-coefficient (D) QSAR / structure-property models for skin & water

Survey of 13 published models that predict the diffusion coefficient D for
small organic molecules in skin compartments (stratum corneum bulk and its
sub-phases, viable epidermis, dermis) and in aqueous solution.

Sister surveys cover Kp permeability and K partition models; this folder
focuses **only on D-predicting models**. Wherever a model fundamentally
predicts P (permeability), the D contribution is back-calculated and noted.

## Comparison table

Sorted by layer, then year.

| Model                              | Year | Layer       | D-form (scaling vs. solute size)                                     | Training n  | Domain                           | DOI |
|------------------------------------|------|-------------|----------------------------------------------------------------------|-------------|----------------------------------|-----|
| Wilke–Chang                        | 1955 | aqueous     | D ~ V_A^(-0.6); D ~ MW^(-0.6) (effective)                            | ~250 pairs  | dilute non-electrolyte           | 10.1002/aic.690010222 |
| Othmer–Thakar                      | 1953 | aqueous     | D ~ V_A^(-0.6) / eta^1.1                                             | ~40         | dilute aqueous                   | 10.1021/ie50519a036 |
| Reddy–Doraiswamy                   | 1967 | aqueous     | D ~ M_B^0.5 / (eta * (V_A V_B)^(1/3)); D ~ MW^(-1/3)                 | ~96 pairs   | binary liquids                   | 10.1021/i160021a012 |
| Hayduk–Laudie                      | 1974 | aqueous     | D = 13.26e-5 * eta^(-1.14) * V_A^(-0.589); D ~ MW^(-0.6)             | 87          | dilute aqueous                   | 10.1002/aic.690200329 |
| Kasting–Smith–Anderson (1992) /    | 1992 | SC          | log10(D/h^2) = a - b*MW; b ~ 0.012-0.018 1/Da; D ~ 10^(-b*MW)        | ~50-80      | human SC, MW 18-400              | book chapter |
| Kasting (finite dose)              | 2001 | SC          | log10(D/h^2) = a - b*MW; b ~ 0.014 1/Da                              | small + 50  | finite-dose human SC             | 10.1002/1520-6017(200102)90:2<202::AID-JPS9>3.0.CO;2-K |
| Mitragotri (SPT)                   | 2002 | SC_lipid    | D_lip ~ exp(-0.46 * r^2);  area-exponential form (r in Angstrom)     | ~60         | hydrophobic, MW < 500            | 10.1002/jps.10031 |
| Mitragotri (4-pathway)             | 2003 | SC_lipid    | D_freeVol ~ exp(-pi (r/r_f)^2) with r_f ~ 4 A; D_lat ~ Saffman-Delbrueck | ~120     | hydrophobic + hydrophilic SC     | 10.1016/S0168-3659(02)00321-8 |
| Anissimov–Roberts (variable-D)     | 2004 | SC          | D(z) functional forms (linear / exp / two-slab); fitted per-compound | n/a (per-cmpd) | depth-profile / desorption fits | 10.1002/jps.10567 |
| Wang–Kasting–Nitsche (D_lip-lat)   | 2007 | SC_lipid    | D_lip-lat = 8.98e-3 * MW^(-2.43) + 2.34e-9                           | ~120 (DB)   | continuous lipid mortar          | 10.1002/jps.20883 |
| Wang–Kasting–Nitsche (k_trans)     | 2007 | SC_lipid    | log10 k_trans = -0.725 - 0.792 * MW^(1/3)  (mass-transfer, cm/s)     | ~120 (DB)   | trans-bilayer hopping            | 10.1002/jps.20883 |
| Lian–Chen–Han (model evaluation)   | 2008 | SC          | benchmark of 7 SC models on 124-compound dataset                     | 124         | SC permeability comparison       | 10.1002/jps.21074 |
| Kretsos–Miller–Zamora–Kasting      | 2008 | dermis      | D_de = D_free(MW) / (1 + B);  log D_free = c0 - c1 * MW^(1/3)        | 26          | mammalian dermis, 4 species      | 10.1016/j.ijpharm.2007.06.020 |
| Chen–Lian–Han                      | 2010 | SC          | combines Mitragotri D_lip + Bungay-Brenner D_cor + tortuous D_pore   | 124         | hydrophilic + hydrophobic SC     | 10.1002/aic.12048 |
| Nitsche–Frasch (binding D_eff)     | 2011 | corneocyte  | D_eff = D_free / (1 + R_B); R_B from keratin binding QSAR            | ~150 (Hansen DB) | corneocyte + binding         | 10.1016/j.ces.2011.01.024 |
| Nitsche–Kasting (corneocyte)       | 2013 | corneocyte  | D_cor = D_aq * (1-phi_f) * H(lambda); Bungay-Brenner polynomial      | 4 (illustr.) | hydrated keratin matrix         | 10.1016/j.bpj.2013.04.009 |
| Nitsche–Kasting (VE)               | 2013 | VE          | D_VE_eff = D_aq / H_epi; H_epi ~ 5-10                                | 4 (illustr.) | viable epidermis                 | 10.1016/j.bpj.2013.04.009 |

## MW-scaling cheat-sheet

The most important number in any of these models is the molecular-weight
exponent (or its analog). They differ wildly:

| Family                          | D ~ MW^? form                | Comment |
|---------------------------------|------------------------------|---------|
| Stokes–Einstein                 | MW^(-1/3)                    | r ~ MW^(1/3); shallow |
| Wilke–Chang / Hayduk–Laudie     | MW^(-0.6)                    | empirical aqueous |
| Reddy–Doraiswamy                | MW^(-1/3)                    | matches Stokes-Einstein |
| Mitragotri 2002 / 2003 lipid    | exp(-c r^2) ~ exp(-c MW^(2/3)) | very steep at large MW |
| Wang–Kasting D_lip-lat          | MW^(-2.43) + floor           | very steep, with non-zero floor |
| Wang–Kasting k_trans            | 10^(-A MW^(1/3)) ~ exp(-c MW^(1/3)) | exponential of cube-root |
| Kasting 1992 / 2001 SC          | 10^(-b MW)  with b ~ 0.014   | exponential in MW (raw) |
| Kretsos D_dermis                | 10^(-c1 MW^(1/3)) / (1 + B)  | exponential cube-root + binding |
| Nitsche-Kasting D_cor / D_VE    | inherit from D_aq            | MW^(-0.6) or MW^(-1/3) via D_aq |

For the **same MW = 200 Da, log K_ow = 2** test compound the predicted
D_SC values across these models span ~5 log decades (1e-7 to 1e-12 cm^2/s)
— the variation comes almost entirely from the MW-exponent assumption,
not from absolute prefactor calibration. Cross-validating with multiple
models gives a useful uncertainty bound for forward-prediction work.

## How to use this folder

Each `*.md` file is a stand-alone model card with frontmatter (`name`,
`predicts`, `layer`, `units`, `descriptors`, `domain`, `status`) and
sections for the citation, equation, training set, performance,
limitations, and references. The frontmatter is intended to be machine-
readable so that a future R helper can iterate over the folder and
build a typed catalog of D predictors that wrap into `skindiff::layer()`'s
D argument.

Status values:
- **classical** — the canonical/historical model in its area (Mitragotri,
  Kasting, Wilke-Chang, Hayduk-Laudie, Othmer-Thakar, Anissimov-Roberts,
  Reddy-Doraiswamy).
- **modern** — current standard for new mechanistic skin models
  (Wang-Kasting, Nitsche-Kasting, Nitsche-Frasch, Kretsos, Chen).
- **benchmark** — comparative evaluation papers, not themselves a
  D model (Lian-Chen-Han 2008).
