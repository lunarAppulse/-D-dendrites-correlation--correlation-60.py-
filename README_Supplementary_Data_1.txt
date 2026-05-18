Supplementary Data 1: Phase-excluded TreeNet dataset

Manuscript:
Interpretable Parameter-Space Mapping of Multiplexed Coding in a Two-Compartment Burst-Propagation Model

Description:
This CSV file contains the phase-excluded dataset used for the TreeNet regression analysis. Each row corresponds to one simulated parameter configuration of the two-compartment burst-propagation model (BPM). The original compiled table contained 1000 simulated configurations. Rows with undefined correlation-based output metrics were excluded before TreeNet fitting, leaving 990 complete observations.

Columns:
sigma: noise amplitude applied to the model input, in pA.
Id0: dendritic input amplitude, in pA.
Is0: somatic input amplitude, in pA.
gs: active coupling to the soma, in model units / pA as used in the simulation code.
cd: active coupling to the dendrite, in model units / pA as used in the simulation code.
Corr_Vs_Vd: population-averaged absolute Pearson correlation between somatic and dendritic membrane potentials.
Corr_ER_Is: Pearson correlation between event rate and somatic input.
Corr_BF_Id: Pearson correlation between burst fraction and dendritic input.
Multiplex_metric: multiplexing metric, defined as the average of Corr_ER_Is and Corr_BF_Id.

Notes:
- The phase-shift parameter was excluded from this TreeNet dataset and was not used as a predictor.
- Undefined output rows were caused by degenerate traces for which one or more Pearson correlations could not be computed.
- The cleaned dataset contains 990 complete rows and 9 columns.
