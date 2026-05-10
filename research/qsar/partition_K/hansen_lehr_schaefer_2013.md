---
name: hansen_lehr_schaefer_2013
predicts: K
phase_pair: SC_lipid/water
units: dimensionless (volume-based)
descriptors: logKow
domain: 16 published K_lip/w measurements + ~75 keratin-binding compounds; log Kow approx -3 to +6
status: modern (review/recommended-parameter paper)
---

# Hansen, Lehr & Schaefer (2013) -- recommended K_lip/water and K_kw correlations

## Citation
Hansen S, Lehr CM, Schaefer UF. Improved input parameters for diffusion models of skin absorption. Adv Drug Deliv Rev. 2013 Feb;65(2):251-264. DOI: 10.1016/j.addr.2012.04.011. PMID: 22626979.

## Equation
The recommended SC lipid/water partition correlation (state-of-the-art at time of publication, fit on 16 carefully measured K_lip/w data points):

    log K_lip/w = 0.092 + 0.67 * log K_ow
    -->  K_lip/w  =  1.23 * K_ow^0.67     (i.e. (alpha, beta) = (0.092, 0.67) in the log-log sense)

so the **K-vs-Kow exponent is 0.67** for the lipid phase. The companion keratin-binding correlation (taken from the Hansen et al. 2011 extended database, 64-compound consolidated set):

    log K_kw = 0.65 + 0.32 * log D       (D = pH-corrected K_ow at pH 7.4)

with slope ~0.32 (keratin-water binding is much less Kow-sensitive than lipid partitioning, in line with the Anderson-Raykar/Nitsche-Wang-Kasting protein exponent of 0.31).

## Training set
- K_lip/w fit: 16 published measurements where the lipid phase was directly extracted/separated; spans log Kow approx -2 to +5.
- Keratin K_kw fit (Hansen 2011 part): 64 compounds in the unified database after curation; log D range -4.67 to +5.70, MW 32-1373 Da.

## Reported performance
- K_lip/w fit: not stated in this 2013 review, but inherits the scatter of the 16 original measurements (typically rms log10 K_lip/w ~ 0.3-0.4).
- K_kw (Hansen 2011): R^2 = 0.78, SE = 0.35.

## Validity / limitations
- Hydrated human SC, lipid phase isolated by extraction; in-tissue value of K_lip/w can differ from the extract-phase value because the bilayer order is lost on extraction.
- Lipid prefactor 1.23 (not 0.43 like Wang/Kasting) reflects the choice of "lipid" reference: extract-phase K vs in-bilayer K. Important: do not blindly mix the Hansen prefactor with Wang/Kasting volume-fraction sums.
- Designed as input to mechanistic skin-absorption transport models; the authors explicitly discuss limitations of single-power-law fits and recommend using both K_lip and K_kw in a two-domain transport simulator.

## Notes
This is the recommended K parameterization for `skindiff`-style multi-layer transport models *as of 2013*: a deliberately curated, separately measured K_lip/w slope, plus a separately measured K_kw slope. The 0.67 exponent sits in the middle of the range across the literature (Cleek-Bunge 0.74, Nitsche-Wang-Kasting 0.81, Frasch-Barbero 0.69, Mitragotri 0.7). The 0.32 keratin-binding slope is now the consensus value across multiple recalibrations.

## References
- Hansen, Lehr, Schaefer 2013 doi:10.1016/j.addr.2012.04.011
- Hansen, Selzer, Schaefer, Kasting 2011 (keratin database) doi:10.1002/jps.22396
