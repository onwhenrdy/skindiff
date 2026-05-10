---
name: mitragotri_2003
predicts: K
phase_pair: SC_lipid/water (and via four-pathway sum, K_SC/water)
units: dimensionless (volume-based, lipid bilayer reference)
descriptors: logKow, molecular radius (Angstrom)
domain: small organic permeants, MW < 500 Da; mostly log Kow -2 to +5; aqueous donor solutions
status: classical (mechanistic)
---

# Mitragotri (2002 / 2003) -- Scaled Particle Theory K_lip/water

## Citation
- Mitragotri S. A theoretical analysis of permeation of small hydrophobic solutes across the stratum corneum based on Scaled Particle Theory. J Pharm Sci. 2002 Mar;91(3):744-752. DOI: 10.1002/jps.10048. PMID: 11920759.
- Mitragotri S. Modeling skin permeability to hydrophilic and hydrophobic solutes based on four permeation pathways. J Control Release. 2003 Jan 9;86(1):69-92. DOI: 10.1016/S0168-3659(02)00321-8. PMID: 12490374.

## Equation
Scaled Particle Theory (SPT) predicts that the lipid-bilayer/water partition coefficient scales as the octanol/water partition coefficient with an exponent slightly less than 1 because both phases are hydrophobic but the lipid bilayer is more ordered. The permeability through the lipid pathway is

    P_lipid  =  5.6e-6 * K_ow^0.7 * exp(-0.46 * r^2)        [cm/s]

with r the solute hard-sphere radius in Angstrom. The K piece is

    K_lip/w  ~  K_ow^0.7

so the **K-vs-Kow exponent in the Mitragotri SPT model is ~0.7**. The 2003 paper extends to four parallel pathways (free-volume diffusion through lipid bilayers, lateral diffusion along bilayers, polar pores, lipid shunts) but the K_lip/w exponent stays ~0.7 across them.

## Training set
The 2003 paper validates against ~100+ compounds from the Flynn 1990 database plus subsequent measurements; the K-side parameterization was inherited from earlier liposome-bilayer partitioning data and SPT theoretical estimates rather than fit from scratch.

## Reported performance
"Best predictor of skin permeability among 7 models compared" in independent benchmarks (Mitragotri et al. 2011, Lian et al. 2008 evaluation). RMSE on log K_p approximately 0.6-0.7 against the standard skin-permeability test sets.

## Validity / limitations
- Solutes < 500 Da, with hard-sphere radius computable from MW or molecular surface.
- Only neutral (uncharged) species; ionizable compounds need fraction-neutral correction at the relevant pH.
- The K_lip/w piece predicts uniformly *higher* K than the Anderson-Raykar / Nitsche-Wang-Kasting recalibrations (different prefactor and slightly different slope), so the absolute lipid-phase K from this model and from the Wang/Kasting model can disagree by a factor of 2-3 even when the Kow exponents are similar.

## Notes
The Mitragotri model is "mechanistic" rather than empirical: K_lip/w is *derived* from SPT estimates of cavity-formation work in an ordered lipid bilayer plus solute-lipid interaction, not regressed against an SC partition database. This makes the K vs Kow slope close to but not exactly 0.7 for any given molecule -- there is a weak size dependence absorbed into the exponent across the ensemble. For `skindiff` purposes, treat 0.7 as the model's K-vs-Kow exponent.

## References
- Mitragotri 2002 doi:10.1002/jps.10048
- Mitragotri 2003 doi:10.1016/S0168-3659(02)00321-8
