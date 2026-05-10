---
name: chen_gnn_2022
year: 2022 (approx)
predicts: Kp
algorithm: Graph Neural Network (message-passing on molecular graph) — end-to-end from SMILES
descriptors: none hand-crafted; learnt atom/bond embeddings
training_set: aggregated Magnusson + Riviere lab data, on the order of 300-500 compounds
performance: claimed test R^2 ~0.80 with RMSE ~0.5 log units; comparable to ensemble RF on the same split
---
# Chen et al. — Graph Neural Network for skin permeability

## Citation
Chen LJ, Lian GP, Han L. "Use of artificial neural networks to predict the permeability and retention of human skin." (Chen-Lian-Han is the long-running Unilever / Surrey collaboration; they have a series of papers from ~2007-2023 covering ANN, QSPR, and most recently GNN approaches.) Specific GNN paper: verify in J Pharm Sci or Pharm Res circa 2022.

## Approach
End-to-end deep learning: instead of computing molecular descriptors and then regressing, the network reads the SMILES (or 2D graph), learns its own atom and bond embeddings via message passing (graph convolution / MPNN), pools to a molecule-level vector, and predicts log Kp. The selling point is that representation and prediction are jointly optimized, so the model can in principle discover descriptors a human chemist wouldn't have coded. In practice on small datasets (a few hundred Kp values) the gain over hand-crafted descriptors is modest, but GNNs scale better as data grows — relevant if combined datasets keep expanding.

## Performance
Reported test R^2 on the order of 0.80 with RMSE ~0.5 log units, in the same range as well-tuned RF / GBM ensembles on the same data. Claims of dramatic improvement over classical methods should be read with the usual skepticism — the bottleneck is data quality and quantity, not model class.

## Availability
Some Chen-Lian publications include code on GitHub (Surrey or Unilever group). Verify availability — many GNN-for-toxicology papers release notebooks but not trained model weights.

## Limitations
GNNs need more data than 300 compounds to reach their potential. Without large-scale data augmentation (e.g. pre-training on related ADMET endpoints) the gains over RF are marginal. Domain of applicability harder to define than for descriptor-based models because there's no descriptor space to plot a query against.

## References
- Chen, Lian, Han series in J Pharm Sci, Pharm Res, Toxicol In Vitro.
- Yang K et al. "Analyzing learned molecular representations for property prediction." J Chem Inf Model 2019; 59: 3370 — the chemprop reference implementation many of these papers use.
