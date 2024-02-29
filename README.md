# Hodgkin EBV MIBI

## Abstract

Classic Hodgkin Lymphoma (cHL) is a tumor composed of rare malignant Reed-Sternberg (RS) cells nested within a tumor microenvironment (TME) consisting of extensive but ineffective immune infiltrates. cHL is associated with Epstein-Barr Virus (EBV) in around 25\% of all cases. The specific contributions of EBV to the pathogenesis of cHL remain largely unknown, in part due to technical challenges in assessing the TME in high detail. Herein, we applied multiplexed ion beam imaging (MIBI) spatial proteomics on 6 EBV-positive and 14 EBV-negative cHL samples. We identify key TME features that distinguish between EBV-positive and EBV-negative cHL, including the predominance of memory CD8 T cells and increased T cell dysfunction as a function of spatial proximity to RS tumor cells. Building upon a larger multi-institutional cohort of 22 EBV-positive and 24 EBV-negative cHL samples, we orthogonally validated our findings through a spatial multi-omics approach, coupling whole transcriptome capture with antibody-defined cell types for tumor and T cell populations within the cHL TME. We delineate contrasting transcriptomic immunological signatures between EBV-positive and EBV-negative cases that differently impact RS cell proliferation, tumor-immune interactions, and mechanisms of T cell dysfunction and dysregulation. Our multi-modal framework enabled a comprehensive dissection of EBV-linked reorganization and immune evasion within the cHL TME, and further highlights the need to elucidate the cellular and molecular factors of virus-associated tumors, with potential for targeted therapeutic strategies.

## Overview of code

The following is an overview of the code in this repository. Please note that since this repository doesn't contain the necessary data that was used for the analysis, if one wishes to reproduce the result, one will need to download the data and change the file path accordingly in the code.

| File Name | Description |
| :---------- | :---------- |
| Figure2.rmd | Produces Figure 2A, 2B, 2C, 2E   |
| Figure3.rmd | Produces Figure 3A, 3B, 3C, 3D, 3E            |
| Figure4.rmd | Produces Figure 4B, 4C           |
| Figure5/ | Includes scripts that produce Figure 5, and related supplementary figures |
| Figure6.rmd | Produces Figure 6A, 6B, 6C, 6D |
| Supplementary.rmd| Produces supplementary figures for Figures 2, 3, and 4|
| multi_test_correction.R | Conducts all the statistical tests in the paper and corrects for multiple comparisons using the Benjamini-Hochberg procedure |
|Batch effect correction and evaluation.R| This script is to adjust and assess data for batch effects based on GeoMX data|

