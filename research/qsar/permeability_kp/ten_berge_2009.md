---
name: ten_berge_2009
predicts: Kp (split into lipid and protein/keratin pathways)
units: cm/h
descriptors: logKow, MW (also vapor pressure and water solubility for the full IH-SkinPerm flux model)
domain: occupational risk assessment; aqueous and vapor donors; MW 18-584, log Kow -3.7 to +5.5
status: classical
---

# ten Berge (2009) — SKINPERM / IH-SkinPerm

## Citation
ten Berge, W. (2009). QSARs for skin permeation of chemicals.
http://home.planet.nl/~wtberge/qsarperm.html (also distributed as IH-SkinPerm
v2.0, May 2017, by AIHA / Tibaldi, ten Berge, Drolet).

Companion: Tibaldi, R., ten Berge, W., Drolet, D. (2014). Dermal absorption
of chemicals: estimation by IH SkinPerm. J. Occup. Environ. Hyg.
11(1):19-31. DOI: 10.1080/15459624.2013.831983.

## Equation
SKINPERM splits Kp into two parallel pathways:

Kp_total = Kp_lipid + Kp_protein   (cm/h, in series with viable epidermis)

Kp_lipid is a Potts-Guy-style nonlinear regression of measured aqueous Kp
against log Kow and MW; Kp_protein (also called Kp_keratin) is a separate
small-residual-flux pathway dominant for hydrophilic small molecules. The
exact regression coefficients live inside the SKINPERM Excel macros and are
not published as a single closed-form equation in the original ten Berge
paper. Practical use is via the Excel tool or its R/Python ports.

For the lipid pathway, ten Berge's regression form (Tibaldi et al. 2014
report logical structure, exact coefficients in the workbook) is:

log Kp_lipid [cm/h] = b0 + b1 * log Kow + b2 * MW + b3 * (log Kow)^2 + b4 * (logKow * MW)

with b0 ~ 4.21, b1 ~ -1.74, b2 ~ -0.060, b3 ~ 0.72, b4 ~ 0.00030 (best-fit
coefficients reported by Tibaldi et al. for aqueous donor; verify against
the workbook before relying).

## Training set
n = 182 measured Kp values for substances with log Kow between -4.49 and
+6.13. Test set of 27 structures used for external validation. Built on the
EDETOX / Vecchia-Bunge / Patel compilation, which subsumes Flynn 1990 and
adds post-1990 measurements.

## Reported performance
Predicted log Kp values mostly within one order of magnitude of measured
values. R^2 ~ 0.78-0.85 reported for various model variants. The
two-pathway split is the main accuracy gain over Potts-Guy.

## Validity / limitations
- Substance domain explicitly stated: MW 18-584 Da; log Kow -3.7 to +5.5.
  Outside this range the model extrapolates but ten Berge cautions against
  it.
- Salts of strong acids/bases: enter -3 as the log Kow.
- Vapor exposures handled separately via a stagnant-air-layer model.
- Maximum dermal absorption is bounded by water solubility * Kp_aq -- the
  tool issues a warning when this bound is approached.
- Includes a viable-epidermis correction analogous to Cleek-Bunge.

## Notes
The IH-SkinPerm Excel tool is the workhorse for occupational dermal-risk
assessment in Europe and (via AIHA) in the US. It packages a Potts-Guy /
Cleek-Bunge / two-pathway / vapor-deposition model into one workflow. For
QSAR-only Kp prediction (without the full kinetic simulation) users should
extract the SKINPERM lipid+protein regression form. The exact coefficients
are not in a peer-reviewed paper -- the only authoritative source is the
ten Berge website / workbook.

## References
- ten Berge, W. (2009). QSARs for skin permeation of chemicals.
  http://home.planet.nl/~wtberge/qsarperm.html
- Tibaldi, R., ten Berge, W., Drolet, D. (2014). Dermal absorption of
  chemicals: estimation by IH SkinPerm. J. Occup. Environ. Hyg. 11(1):19-31.
  DOI: 10.1080/15459624.2013.831983.
- IH SkinPerm v2.0 Reference Manual (May 2017): AIHA Exposure Assessment
  Strategies Committee.
