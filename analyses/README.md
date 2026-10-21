# Analyses

This directory contains the scripts used for downstream analyses and figure generation. The tables below indicate the correspondence between scripts and figures in the manuscript.

## Main analyses

| Script | Figure(s) |
|---|---|
| `aba_mutants_analysis.R` | Fig. 4d–e; Extended Data Fig. 11a |
| `dapseq_targets.R` | Fig. 4b (prediction of targets for regulon analysis) |
| `drought_composition_analysis.R` | Fig. 5e; Extended Data Fig. 12d |
| `drought_embedding.R` | Fig. 5d |
| `drought_expression_levels.R` | Extended Data Fig. 12e |
| `drought_fluorecsence_per_cell_ecdf.R` | Fig. 5h |
| `drought_umi_adjusted.R` | Extended Data Fig. 12f–g |
| `dual_marker_quants.R` | Fig. 3e, g, h |
| `ecotype_deg.R` | Extended Data Fig. 8e–f |
| `ecotype_ERS_prop.R` | Extended Data Fig. 8d |
| `ecotype_module_scores.R` | Extended Data Fig. 8b–c |
| `hdWGCNA_module_detection.R` | Fig. 2b–c |
| `hdWGNCA_regulons_2.R` | Fig. 4b |
| `MS_concentration_analysis.R` | Fig. 5a–c |
| `PAGA_trajectories.R` | Fig. 2e–g |
| `plot_ERS_features.R` | Fig. 2d; Fig. 3a |
| `SCP_cluster_markers_plot.R` | Fig. 1c |
| `scr3_analysis.R` | Extended Data Fig. 11b–d |

## Quality control (`QC/`)

The `QC/` directory contains scripts used for various quality control analyses, including comparisons with published datasets, assessment of cell quality and evaluation of developmental trajectories.

| Script | Figure(s) |
|---|---|
| `cell_cycle.R` | Extended Data Fig. 3a–c |
| `cellranger_mapping_wendrich.sh` | Data processing for Extended Data Fig. 6 |
| `comparison_farmer.R` | Extended Data Fig. 5a–f |
| `comparison_shahan.R` | Extended Data Fig. 7a–f |
| `comparison_wendrich.R` | Extended Data Fig. 6a–g |
| `correlation_pseudotime_shahan_marker_scores.R` | Extended Data Fig. 2a |
| `downsampling.R` | Extended Data Fig. 4d–g |
| `expression_levels_ERS.R` | Extended Data Fig. 4a–c |
| `farmer_id_projection_stability.R` | Extended Data Fig. 5g–k |
| `module_genes_de_protoplasting.R` | Extended Data Fig. 3e |
| `monocle_analysis.R` | Extended Data Fig. 2b–e |
| `percent_mito_chloro_genes.R` | Extended Data Fig. 3i–j |
| `scvelo.R` | Extended Data Fig. 2f–g |
| `shahan_id_projection_stability.R` | Extended Data Fig. 7g–k |
| `through_agar_analysis.R` | Extended Data Fig. 3f–h |
| `wendrich_id_projection_stability.R` | Extended Data Fig. 6h–l |

## snATAC-seq analysis (`snATAC_Farmer/`)

This directory contains scripts for the processing and analysis of the single-nucleus ATAC-seq dataset published by Farmer et al. (2021), as presented in Extended Data Fig. 9.
