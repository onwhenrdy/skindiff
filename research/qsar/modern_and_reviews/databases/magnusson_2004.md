---
name: magnusson_2004
type: database
url: not directly hosted; as supplementary table in Magnusson 2004 J Invest Dermatol paper
license: free (publication-derived)
---
# Magnusson 2004 — extended Kp dataset

## Description
Brett Magnusson, Hadgraft, and colleagues extended the Flynn 1990 dataset to ~278 compounds with measured human or animal Kp values, and re-fit a maximum-flux QSAR (J_max instead of just Kp). This is the second-most-cited skin permeation dataset and the standard "extended training set" for modern QSAR / ML papers.

## Coverage
- ~278 compounds in the J_max set; corresponding Kp values for the subset where they're available.
- Mix of human cadaver skin and porcine / hairless mouse skin (with skin-source corrections).
- Adds many drug-like molecules to the Flynn small-organic core: NSAIDs, beta-blockers, antifungals, peptides.
- Vehicle still mostly aqueous; some entries from non-aqueous vehicles flagged.
- Kp range similar to Flynn (~5 log units); J_max adds molecular-weight-dependent solubility ceiling info.

## Access
Supplementary material of Magnusson BM, Anissimov YG, Cross SE, Roberts MS. "Molecular size as the main determinant of solute maximum flux across the skin." J Invest Dermatol 2004; 122(5): 993-999.

DOI: 10.1111/j.0022-202X.2004.22413.x.

The compound list appears in subsequent QSAR papers (Sun-Moss-Davey series, Baba 2015, Tsakovska 2017) often with curation flags.

## Citation
Magnusson BM, Anissimov YG, Cross SE, Roberts MS. "Molecular size as the main determinant of solute maximum flux across the skin." J Invest Dermatol 2004; 122: 993-999.

## Notes
- The Magnusson J_max model: log J_max ~ - 3.90 - 0.0190 * MW (in J_max in mg/cm^2/h). Simple but surprisingly hard to beat.
- Together with Flynn, this is the standard "skin Kp benchmark" - any new QSAR's headline number is its R^2 on (some subset of) Flynn + Magnusson.
- Same fundamental data quality issues as Flynn: heterogeneous skin sources, lab-to-lab variation, and the noise floor of in-vitro permeation experiments (~0.3-0.5 log units for the same compound across labs).
