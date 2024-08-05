# The Spontaneous Neoantigen-Specific CD4+ T Cell Response to a Growing Tumor is Functionally and Phenotypically Diverse

This repository contains all the necessary resources and scripts used in our study, **"The Spontaneous Neoantigen-Specific CD4+ T Cell Response to a Growing Tumor is Functionally and Phenotypically Diverse."** Below are the details of the tools and data involved.

## Software and Modules

This project utilizes the following software and versions:

- **R**: v3.6.1
- **Cell Ranger**: v3.1.0
- **Seurat**: v3.1.5

## Raw Data

The raw and processed single-cell RNA-seq and TCR-seq files are available via the GEO accession number.

## Data Pre-processing

### Single-cell Analysis

- **Demultiplexing and Mapping**: Utilize our in-house pipeline along with Cell Ranger for demultiplexing and mapping.
- **Quality Control**: Employ our in-house pipeline for single-cell quality control measures.
- **Clustering**: Generate clusters using our in-house pipeline integrated with Seurat.
- **VDJ Libraries Aggregation**: Aggregate VDJ libraries using our in-house [pipeline](https://github.com/vijaybioinfo/VDJ_aggr).

## Citation

If you utilize this repository for your research, please cite our study:

**The Spontaneous Neoantigen-Specific CD4+ T Cell Response to a Growing Tumor is Functionally and Phenotypically Diverse**

## Contact

For further inquiries or assistance, please reach out to:

- **Manuel Azaid Ordaz Arias**: [moarias@lji.org](mailto:moarias@lji.org)
- **Ryan Griswold**: [rgriswold@lji.org](mailto:rgriswold@lji.org)
- **Vijayanand Pandurangan**: [vijay@lji.org](mailto:vijay@lji.org)
