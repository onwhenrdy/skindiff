---
name: chemaxon_calculator_plugins
type: software package (commercial)
url: https://chemaxon.com/products/calculators-and-predictors
license: commercial (free for academic use under their academic license)
---
# ChemAxon Calculator Plugins

## Description
ChemAxon's Calculator Plugins are a Java-based suite of property predictors integrated into their JChem / Marvin chemistry framework. They cover ~40 physicochemical properties with strong descriptor coverage: logP / logD profile, pKa, polar surface area, hydrogen bond donors/acceptors, refractivity. Skin permeability is offered as a derived property using a Potts-Guy-style equation on the calculated descriptors.

Stronger areas than skin Kp specifically: logP / logD, pKa (the cxcalc pKa predictor is widely regarded as one of the best free / low-cost pKa tools).

## Coverage
- ~40 physicochemical properties.
- Skin Kp: derived from logP and MW, no separate training data — equivalent to applying Potts-Guy to ChemAxon-calculated logP.
- High throughput, batch-friendly: tens of thousands of compounds per second on common hardware.

## Access
- Commercial; pricing varies by deployment model. Academic licensing is generous (free for many academic uses).
- API: Java library, command-line (cxcalc), KNIME nodes, Pipeline Pilot integration, R/Python wrappers.
- Web service / REST API also available with a hosted license.

## Citation
ChemAxon Ltd, Budapest. *Marvin / JChem Documentation*. https://chemaxon.com/

No primary academic citation for the proprietary methods; the underlying logP method is documented as "VG / Klopman fragmental approach with corrections."

## Notes
- ChemAxon's logP / pKa predictors are commonly used as descriptor inputs for downstream skin Kp QSARs (Potts-Guy or ML).
- Built-in skin permeability is a thin wrapper around Potts-Guy and adds little vs. running Potts-Guy directly with any logP source.
- API quality and language coverage make ChemAxon a common choice for industrial pipelines that need batch property calculation.
