---
name: oecd_qsar_toolbox
type: toolbox / software package
url: https://qsartoolbox.org/
license: free (OECD-supported public tool)
---
# OECD QSAR Toolbox

## Description
The OECD QSAR Toolbox is a free desktop application developed jointly by the OECD and the European Chemicals Agency (ECHA), with software development by Bulgarian Laboratory of Mathematical Chemistry (LMC) at the University of Burgas. It is a hybrid toolbox combining: (a) QSAR-grouping (read-across) workflows, (b) a built-in database of measured property data including dermal absorption, and (c) several mechanistic and statistical models. Used by industry and regulators globally for filling data gaps in chemical safety dossiers.

## Coverage
- Hundreds of thousands of compounds with at least one measured property.
- Dermal-absorption-relevant properties: in-vitro Kp (Flynn / Magnusson / additions), some in-vivo dermal absorption fractions, octanol-water logP, water solubility, MW.
- Profilers for skin sensitization (the toolbox's strongest dermal area), skin irritation, and metabolic activation.
- Built-in skin Kp QSARs include Potts-Guy and several variants.

## Access
Free download from https://qsartoolbox.org/ after registration. Windows desktop application; some Linux support via wine. Periodic version updates (current major version 4.x as of recent years).

Includes API for batch processing via command-line interface. Limited scripting / programmatic access compared to a true library.

## Citation
Diderich R. "OECD QSAR Toolbox." In: Cronin MTD, Madden JC, Enoch SJ, Roberts DW, eds. *Chemical Toxicity Prediction: Category Formation and Read-Across*. RSC 2013.

Plus annual OECD QSAR Toolbox release notes, available from the website.

## Notes
- Primary strength is the read-across / category-grouping workflow for regulatory dossier preparation, not as a Kp prediction engine.
- Skin permeability models are present but limited; for serious Kp prediction, dedicated tools (Tsakovska consensus, COSMOtherm-Skin, ML papers) are better.
- The dermal data subset of the database is one of the better-curated public collections of in-vitro Kp values, useful as a starting point for any new training set.
