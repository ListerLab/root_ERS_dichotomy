# An environmentally regulated dichotomy in cell lineages of the root

This repository contains the code used for the processing and analysis of the single-cell RNA-sequencing datasets presented in our study investigating an environmentally regulated cell state (ERS) across *Arabidopsis thaliana* root cell lineages.

**Oliva M. et al. (2026). An environmentally regulated dichotomy in cell lineages of the root.**

[Link to publication]


## Repository structure

### `data/`

Scripts used for processing sequencing datasets and generating the objects used for downstream analyses.

- `pre-processing_per_library/` – pre-processing and quality control of individual scRNA-seq libraries.
- `integration/` – integration of libraries for the different datasets and experimental comparisons used in the study.

Further information is provided in [`data/README.md`](data/README.md).

### `analyses/`

Scripts used for downstream analyses presented in the manuscript.

A correspondence between analysis scripts and manuscript figures is provided in [`analyses/README.md`](analyses/README.md).


## Data availability

Raw FASTQ files generated in this study have been deposited in the European Nucleotide Archive (ENA) under Project accession [PRJEB100815](https://www.ebi.ac.uk/ena/browser/view/PRJEB100815).

Processed Seurat objects used for downstream analyses have been deposited on Zenodo under the following DOIs:

- [10.5281/zenodo.17411120](https://doi.org/10.5281/zenodo.17411120)
- [10.5281/zenodo.17411065](https://doi.org/10.5281/zenodo.17411065)

Previously published datasets reanalysed in this study are available from their original repositories. Details and accession information are provided in the Methods of the associated manuscript.

## Software

Software and package versions used for data processing and analysis are described in the Methods of the associated manuscript.


## Contact

For questions regarding the code or analyses, please contact:

**Marina Oliva**  
marina.oliva@uwa.edu.au
