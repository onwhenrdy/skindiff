---
name: admetlab
type: web server
url: https://admetmesh.scbdd.com/ (ADMETlab 2.0); https://admetlab3.scbdd.com/ (3.0, current)
license: free for academic use
---
# ADMETlab 2.0 / 3.0

## Description
ADMETlab is a free public web server for ADMET property prediction developed by the Tang lab at Central South University (Changsha, China). It predicts ~90 endpoints including absorption, distribution, metabolism, excretion, and toxicity properties. Skin permeability (logKp from Potts-Guy-style descriptors) is one of the absorption endpoints. Among the most user-friendly and broadly used free ADMET web servers globally.

## Coverage
- ~90 ADMET endpoints in v3.0 (2024).
- Skin Kp: a small ML model (RF / SVM-class) trained on Flynn-Magnusson data.
- High throughput web interface; batch upload supported.

## Access
- Free web interface at https://admetlab3.scbdd.com/ (current generation) or https://admetmesh.scbdd.com/ (2.0, still up).
- Accepts SMILES (single or batch).
- API access via documented REST endpoints; some scripting support.
- No download of trained models.

## Citation
- Xiong G, Wu Z, Yi J, Fu L, Yang Z, Hsieh C, Yin M, Zeng X, Wu C, Lu A, Chen X, Hou T, Cao D. "ADMETlab 2.0: an integrated online platform for accurate and comprehensive predictions of ADMET properties." Nucleic Acids Research 2021; 49(W1): W5-W14. DOI: 10.1093/nar/gkab255.
- v3.0 release paper: Fu L et al. Nucleic Acids Research 2024; 52(W1): W422-W431.

## Notes
- Excellent for fast ADMET screening including a quick skin Kp number for any chemical.
- The Kp model is one of many endpoints and is not the strongest in the suite — for serious skin Kp work, supplement with dedicated tools.
- Server hosted in China; latency from elsewhere is variable. Availability and policy may change.
- Comparable free web servers: ADMET-AI (Stanford / Greg Landrum, https://admet.ai.greenstonebio.com/), pkCSM (https://biosig.lab.uq.edu.au/pkcsm/), Swiss-ADME (https://www.swissadme.ch/) — though the last has no skin Kp endpoint.
