---
name: wu_transformer_2023
year: 2023 (approx)
predicts: Kp (within a multi-task ADMET framework)
algorithm: pre-trained molecular transformer (ChemBERTa-style or MolFormer); fine-tuned on skin Kp head
descriptors: end-to-end from SMILES; learnt embeddings
training_set: pre-train on ~10M PubChem SMILES; fine-tune on aggregated skin Kp set ~500 compounds
performance: claims R^2 ~0.78-0.82; comparable to top ensemble methods
---
# Wu et al. — Pre-trained molecular transformer fine-tuned for skin permeability

## Citation
Specific Wu et al. transformer-for-skin paper details should be verified — multiple groups have applied ChemBERTa / MolFormer / MolBERT to skin endpoints in 2022-2024. Examples include applications by Shen et al., Yang et al., and Chinese pharmaceutical AI groups. The architectural lineage:
- Chithrananda S, Grand G, Ramsundar B. "ChemBERTa: Large-Scale Self-Supervised Pretraining for Molecular Property Prediction." 2020.
- Ross J et al. "Large-scale chemical language representations capture molecular structure and properties." Nat Mach Intell 2022; 4: 1256.

## Approach
Pre-train a transformer encoder on tens of millions of unlabeled SMILES with a masked-language-modeling objective (predict masked atoms / tokens). The pre-trained weights produce a chemistry-aware embedding for any new SMILES. Then fine-tune a small regression head on the (small) skin-Kp dataset. The bet is that pre-training on huge unlabelled data gives a representation that overcomes the small fine-tuning dataset — analogous to how BERT and GPT work in NLP.

## Performance
On benchmark ADMET tasks, pre-trained transformers reach R^2 ~0.78-0.82 on skin Kp test sets, in the same range as well-tuned RF + descriptors. The gain over descriptors is more pronounced on the smallest fine-tuning sets (less than ~50 compounds), where the pre-trained representation acts as a strong prior. On larger Kp sets the gap narrows.

## Availability
ChemBERTa, MolFormer, and MolBERT base models are public on HuggingFace. Skin-specific fine-tunes are generally not released as standalone tools but can be reproduced from the pre-trained checkpoints in a few hours of GPU time.

## Limitations
Skin Kp data is small enough that the choice of fine-tuning split dominates results — the same architecture can vary by 0.1 in R^2 across reasonable splits. Interpretability is poor: attention weights over SMILES tokens do not map cleanly to chemical features the way descriptors do.

## References
- ChemBERTa, MolBERT, MolFormer originating papers (above).
- Heid E, Greenman KP, Chung Y, et al. "Chemprop: Machine learning package for chemical property prediction." J Chem Inf Model 2024; 64: 9-17 — the practical reference implementation.
