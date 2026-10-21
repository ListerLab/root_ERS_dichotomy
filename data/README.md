# Data processing

This directory contains the scripts used for pre-processing and integration of the scRNA-seq datasets.

## Pre-processing (`pre-processing/`)

The pre-processing procedures are described in detail in the manuscript Methods, illustrated in Supplementary Fig. 1, and summarised for each library in Supplementary Table 2.

The following workflow was applied independently to each scRNA-seq library:

```text
┌──────────────────────────────┐
│     Cell Ranger mapping      │
└──────────────────────────────┘
               ↓
┌──────────────────────────────┐
│     Low-quality cell         │
│         filtering            │
└──────────────────────────────┘
               ↓
┌──────────────────────────────┐
│   Doublet identification     │
│         and removal          │
└──────────────────────────────┘
               ↓
┌──────────────────────────────┐
│     SNP-based genotyping     │
│   (mixed-ecotype libraries)  │
└──────────────────────────────┘
```

Scripts applicable to all libraries (mapping, doublet detection, and SNP counting) are located directly in `pre-processing/`. Scripts for low-quality cell filtering and SNP-based clustering, which require library-specific adjustments, are organised in the corresponding library subdirectories.

*Note: Libraries RL970, RL1399, RL1463, RL1464, RL1465, and RL4460 contained mixtures of two ecotypes, one of which was not included in the study. To avoid ambiguity about which cells belong to the published dataset, raw FASTQ files for these libraries are not provided. Instead, the count matrices of the cells included in the study are available in the processed Seurat objects deposited on Zenodo.*

## Integration (`integration/`)

This directory contains the scripts used to integrate biological replicates of each experimental condition and annotate the resulting datasets using the Col-0 reference atlas.

Replicates were integrated using the RPCA integration method implemented in Seurat v5, followed by cell identity annotation through reference-based label transfer.
