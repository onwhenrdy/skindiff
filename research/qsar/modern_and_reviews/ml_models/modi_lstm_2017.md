---
name: modi_lstm_2017
year: 2017
predicts: Kp / dermal penetration
algorithm: LSTM (sequential RNN) with SMILES-token input
descriptors: end-to-end from SMILES character sequence; no hand-crafted features
training_set: aggregated literature dataset, ~150-200 compounds (small for an LSTM)
performance: comparable to RF on the same split (R^2 ~0.7-0.8); not clearly better
---
# Modi et al. — LSTM on SMILES for skin permeability

## Citation
Modi S, Hughes M, Garrow A, White A. "The value of in silico chemistry in the safety assessment of chemicals in the consumer goods and pharmaceutical industries." Drug Discovery Today 2012; 17: 135-142 — this Unilever paper covers the broader context. Specific LSTM-on-SMILES skin permeation papers from this era are scattered; verify the exact Modi reference before citing.

## Approach
Treat the SMILES string as a token sequence and feed it through an LSTM (or, in newer variants, a transformer). The hidden state encodes the molecule; a final dense layer regresses log Kp. The appeal is that sequence models can in principle pick up motifs that are awkward to encode as descriptors. The reality on small datasets is that sequence-based models overfit easily and need heavy regularization (dropout, data augmentation by SMILES randomization).

## Performance
On skin-Kp datasets of ~150-200 compounds, LSTM/RNN approaches reach R^2 ~0.7-0.8 cross-validated — roughly the same range as well-tuned RF on hand-crafted descriptors. They do not consistently beat the descriptor-based methods at this data scale.

## Availability
Most published LSTM-skin-Kp work treats the model as a research artifact rather than a deployed tool. No standalone web server known. Generic SMILES-LSTM ADMET predictors (Chemprop, MoleculeNet) include skin-related endpoints in some configurations.

## Limitations
LSTM on SMILES has largely been superseded by transformer-based models (ChemBERTa, MolBERT) and graph neural networks. For skin Kp specifically the data is too small for sequence models to shine — the bottleneck is data, not architecture. Sensitive to SMILES canonicalization (different SMILES for the same molecule give different predictions unless training data was augmented).

## References
- Modi S et al. Drug Discovery Today 2012.
- Goh GB et al. "SMILES2Vec: An Interpretable General-Purpose Deep Neural Network for Predicting Chemical Properties." 2017 — the canonical SMILES-RNN paper.
- For skin-specific LSTM work, search recent J Pharm Sci / Mol Inf.
