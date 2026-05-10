# QSAR / structure-property models for skin K (partition coefficient)

This folder surveys QSAR / structure-property models that predict the **partition coefficient K** for small organic molecules between skin compartments and a reference phase (water or vehicle). K is the equilibrium ratio of the solute concentration in the skin compartment over the solute concentration in the reference (water) phase: K_SC/water = c_SC / c_water at equilibrium.

`skindiff` takes K as a direct per-layer input. The user may want to set K from a QSAR rather than measure it; this folder summarises the published QSARs available.

## Headline finding: the K-vs-Kow exponent

The most important number across all these models is the **exponent on K_ow** in K_layer/water = (prefactor) * K_ow^(exponent). Across the literature this exponent ranges from **~0.3 (protein/keratin phase)** to **~0.8 (lipid phase)**, with typical *whole-tissue* SC values of 0.67-0.81. Whole-K-vs-Kow slopes:

| Source | Phase | K-vs-Kow exponent | Prefactor |
| --- | --- | --- | --- |
| Potts-Guy 1992 (implicit) | SC/water | 0.71 | (implicit, from K_p) |
| Cleek-Bunge 1993 | SC/water | 0.74 | 1.0 |
| Anderson-Raykar 1988 | SC lipid | ~0.7 (orig); 0.81 (recalibrated) | varies |
| Anderson-Raykar 1988 | SC protein | ~0.3 | ~4.2 |
| Frasch-Barbero 2003 | SC lipid | 0.69 | 1.0 |
| Frasch-Barbero 2003 | SC protein | 0.31 | 4.2 |
| Mitragotri 2002/2003 | SC lipid | 0.7 | (SPT-derived) |
| Nitsche-Wang-Kasting 2006 | SC lipid | 0.81 | 0.43 |
| Nitsche-Wang-Kasting 2006 | SC protein | 0.31 | 4.2 |
| Wang-Kasting-Nitsche 2007 | SC lipid + corneocyte (multiphase) | 0.81 / 0.31 | 0.43 / 4.2 |
| Hansen 2011 (keratin) | keratin/water | 0.32 | (intercept 0.65 in log-log) |
| Hansen 2013 | SC lipid | 0.67 | 1.23 (intercept 0.092) |
| Kretsos-Kasting 2008 | dermis/water | ~0 (small lipid bump only at high Kow) | K ~ 0.7-0.9 baseline |
| Nitsche-Kasting 2013 | viable epidermis/water | ~0 (small lipid bump only at high Kow) | K ~ 1 baseline |
| COSMOmic / MD 2023 | SC lipid | 0.74 (best fit through experimental data) | (in-silico, no closed form) |

Practical takeaways:
- **Lipid phase**: K-vs-Kow exponent in the range 0.67-0.81 across all credible literature; consensus ~0.7-0.75.
- **Protein/keratin phase**: K-vs-Kow exponent ~0.31-0.32 across all credible literature.
- **Dermis and viable epidermis**: K is close to 1 with a small Kow correction at high lipophilicity.
- **Whole-tissue SC**: a single power-law fit (Cleek-Bunge 0.74) is convenient but loses the curvature; the two-domain decomposition (Anderson-Raykar / Wang-Kasting / Hansen) is mechanistically correct.

## Comparison table (grouped by phase pair)

### SC/water (whole tissue, single power law)

| Model | Year | Equation form | Training n | Domain | DOI |
| --- | --- | --- | --- | --- | --- |
| Potts-Guy (implicit K) | 1992 | K_SC/w ~ K_ow^0.71 | ~93 | logKow [-3,+6], MW < 750 Da | 10.1023/A:1015810312465 |
| Cleek-Bunge | 1993 | K_SC/w = K_ow^0.74 | ~93 (Flynn 1990) | logKow [-3,+6], MW < 500 Da | 10.1023/A:1018981515480 |

### SC/water (whole tissue, two-domain decomposition)

| Model | Year | Equation form | Training n | Domain | DOI |
| --- | --- | --- | --- | --- | --- |
| Anderson-Raykar | 1988-89 | K_SC/w = phi_lip*K_lip/w + phi_pro*PC_pro/w + phi_w | hydrocortisone esters + later pooled | logKow [-1,+5] | 10.1023/A:1015956705293 / 10.1023/A:1015989929342 |
| Nitsche-Wang-Kasting | 2006 | K_lip/w = 0.43*K_ow^0.81; PC_pro/w = 4.2*K_ow^0.31 | 72 K_SC/w points | logKow ~8 decades | 10.1002/jps.20549 |

### SC lipid/water

| Model | Year | Equation form | Training n | Domain | DOI |
| --- | --- | --- | --- | --- | --- |
| Frasch-Barbero | 2003 | K_lip/w = K_ow^0.69 | inherited from Anderson-Raykar | hydrated SC | 10.1002/jps.10466 |
| Mitragotri (SPT) | 2002, 2003 | K_lip/w ~ K_ow^0.7 | ~100 (Flynn) | MW < 500 Da, neutral | 10.1002/jps.10048 / 10.1016/S0168-3659(02)00321-8 |
| Wang-Kasting-Nitsche | 2007 | K_lip/w = 0.43*K_ow^0.81 (from NWK 2006) | ~150 K_p set | logKow [-3,+6] | 10.1002/jps.20883 |
| Hansen-Lehr-Schaefer | 2013 | log K_lip/w = 0.092 + 0.67*log K_ow | 16 measured K_lip/w | hydrated SC | 10.1016/j.addr.2012.04.011 |
| Piasentin (COSMOmic/MD) | 2023 | beta = 0.74 best-fit; structure-based de novo | 16-20 | MW < 500 Da | 10.1021/acs.jpcb.2c08566 |

### Corneocyte / keratin / water

| Model | Year | Equation form | Training n | Domain | DOI |
| --- | --- | --- | --- | --- | --- |
| Anderson-Raykar | 1988 | PC_pro/w ~ 4.2*K_ow^0.31 (recalibrated) | hydrocortisone esters + pooled | logKow [-1,+5] | 10.1023/A:1015989929342 |
| Hansen 2011 (extended keratin) | 2011 | log K_kw = 0.65 + 0.32*log D | 64 compounds | logD [-4.7,+5.7], MW 32-1373 | 10.1002/jps.22396 |

### Dermis / water

| Model | Year | Equation form | Training n | Domain | DOI |
| --- | --- | --- | --- | --- | --- |
| Kretsos-Kasting | 2008 | K_de ~ 0.7-0.9 + small phi_lip*K_lip/w correction; explicit fu term | 26 compounds, 4 species | MW 18-476 Da | 10.1016/j.ijpharm.2007.06.020 |

### Viable epidermis / water

| Model | Year | Equation form | Training n | Domain | DOI |
| --- | --- | --- | --- | --- | --- |
| Nitsche-Kasting | 2013 | K_VE/w = 0.879*K_cyt/w + 0.001*K_lip/w + 0.120*K_ext/w | inherited | hydrated VE | 10.1016/j.bpj.2013.03.056 |

### Octanol / water (input to all of the above)

| Predictor | Equation form | Training n | Domain | Reference |
| --- | --- | --- | --- | --- |
| KOWWIN (EPI Suite) | atom/fragment additive | ~13 000 | broad organic | EPA EPI Suite (free) |
| ALOGPS | NN ensemble | ~12 000 | broad organic | doi:10.1021/ci049854h |
| ClogP | fragment + corrections | ~10 000 | drug-like | doi:10.1021/cr00020a002 |
| miLogP | group contribution | ~12 000 | drug-like | molinspiration.com |
| XLogP3 | atom-additive + knowledge | ~8 000 | drug-like; in PubChem | doi:10.1021/ci700257y |

## File index

- [anderson_raykar_1989.md](anderson_raykar_1989.md) -- two-domain decomposition foundation
- [cleek_bunge_1993.md](cleek_bunge_1993.md) -- single-power-law SC/water
- [cosmomic_md_2023.md](cosmomic_md_2023.md) -- in-silico K_lip/w
- [frasch_barbero_2003.md](frasch_barbero_2003.md) -- two-phase K for FEM SC models
- [hansen_2011_keratin.md](hansen_2011_keratin.md) -- extended keratin-binding QSAR
- [hansen_lehr_schaefer_2013.md](hansen_lehr_schaefer_2013.md) -- recommended K_lip/w + K_kw set
- [kretsos_kasting_2008.md](kretsos_kasting_2008.md) -- dermis K
- [logKow_predictors_overview.md](logKow_predictors_overview.md) -- KOWWIN / ALOGPS / ClogP / miLogP / XLogP3 brief
- [mitragotri_2003.md](mitragotri_2003.md) -- Scaled Particle Theory K_lip/w
- [nitsche_kasting_2013_VE.md](nitsche_kasting_2013_VE.md) -- viable epidermis K
- [nitsche_wang_kasting_2006.md](nitsche_wang_kasting_2006.md) -- recalibrated two-phase SC/water (benchmark)
- [potts_guy_1992_implied_K.md](potts_guy_1992_implied_K.md) -- implicit K from Potts-Guy K_p
- [wang_kasting_nitsche_2007.md](wang_kasting_nitsche_2007.md) -- microscopic multiphase SC, K parameterisation

## Recommended defaults for `skindiff`

- **SC as a single homogeneous slab**, just want a K from log Kow: use **Cleek-Bunge 1993** (K_SC/w = K_ow^0.74). Easy, defensible, regulatory-standard.
- **SC split into separate lipid and corneocyte sub-layers**: use **Wang-Kasting-Nitsche 2007** (K_lip/w = 0.43 * K_ow^0.81 and K_cor/w from the additive sum with PC_pro/w = 4.2 * K_ow^0.31), or equivalently **Hansen-Lehr-Schaefer 2013** (K_lip/w = K_ow^0.67 with prefactor 1.23). Both are credible; they differ mainly in the lipid prefactor.
- **Dermis layer**: K_de ~ 0.7-0.9, treat as approximately constant unless log Kow > 3, then add the Kretsos-Kasting 2008 small lipid correction.
- **Viable epidermis layer**: K_VE ~ 1, same caveat as dermis.

When fitting `skindiff` to permeation/penetration data, use these QSAR estimates as priors / starting values for K, not as final values. The spread between QSARs (factor 2-3 in K for highly lipophilic compounds) is a useful prior width.
