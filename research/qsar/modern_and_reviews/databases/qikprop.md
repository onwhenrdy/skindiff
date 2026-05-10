---
name: qikprop
type: software package (commercial)
url: https://www.schrodinger.com/products/qikprop
license: commercial (Schrödinger; academic licensing available)
---
# Schrödinger QikProp

## Description
QikProp is Schrödinger's ADMET property prediction module, integrated into the Maestro molecular modelling environment. It predicts ~50 physicochemical and pharmacokinetic properties from 3D structure, including a "Kp" estimate (skin permeability, predicted log Kp from a Schrödinger-maintained QSAR), QPlogPo/w (logP), QPlogS (water solubility), human oral absorption, and Lipinski/Veber rule compliance.

The skin Kp predictor is one of QikProp's lesser-known endpoints; it uses a small published QSAR similar in spirit to Potts-Guy with proprietary descriptor weighting.

## Coverage
- ~50 ADMET endpoints in total.
- Skin Kp range: typical drug-like molecules; out-of-domain warnings for very large or very polar compounds.
- Throughput: thousands of compounds per minute on a desktop machine.

## Access
- Commercial license through Schrödinger; bundled with the standard Maestro suite.
- Academic licensing available at significantly reduced cost (free for some institutions; site licenses common at universities with structural-biology groups).
- API access via the Schrödinger Python API (Maestro integration), command-line, and KNIME / Pipeline Pilot nodes.

## Citation
Jorgensen WL, Duffy EM. "Prediction of drug solubility from structure." Adv Drug Deliv Rev 2002; 54: 355-366.

For QikProp documentation: Schrödinger LLC, *QikProp 7.x User Manual*.

## Notes
- General-purpose ADMET tool; skin Kp is a peripheral endpoint and not the strongest in the suite.
- Training set / methodology for the skin Kp predictor is not fully transparent (closed source).
- For users already in the Schrödinger ecosystem, fine for screening-level Kp triage. For dedicated dermal work, mechanistic tools (skindiff + a curated Kp QSAR) are more transparent.
- Comparable commercial competitors with skin permeability endpoints: ChemAxon Calculator Plugins, ACD/Labs Percepta, Simulations Plus ADMET Predictor.
