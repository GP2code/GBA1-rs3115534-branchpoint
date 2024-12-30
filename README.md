# African Ancestry Neurodegeneration Risk Variant Disrupts an Intronic Branchpoint in *GBA1*

`GP2 ❤️ Open Science 😍`

[![DOI](https://zenodo.org/badge/741747163.svg)](https://zenodo.org/doi/10.5281/zenodo.10484208)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Updated:** December 2024 (Added link to published manuscript)

## Summary
The following repository encompasses all scripts developed and used in the manuscript titled _**"African ancestry neurodegeneration risk variant disrupts an intronic branchpoint in GBA1"**_. This study explores the disease risk mechanism for the non-coding _GBA1_ rs3115534 variant. 

* Manuscript link: https://www.nature.com/articles/s41594-024-01423-2

This is a follow-up to the first largest genome-wide assessment of Parkinson’s disease in the African and African admixed populations, published here:

> Rizig M*, Bandrés-Ciga S*, Makarious MB*, ..., Singleton AB, Okubadejo N, on behalf of the Global Parkinson’s Genetics Program: _**“Genome-wide Association Identifies Novel Etiological Insights Associated with Parkinson’s Disease in African and African Admixed Populations”**_ *Lancet Neurology (2023)*; https://doi.org/10.1016/S1474-4422(23)00283-1 

### Helpful Links 
- [GP2 Website](https://gp2.org/)
    - [GP2 Cohort Dashboard](https://gp2.org/cohort-dashboard-advanced/)
- [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
    - [Other GP2 Manuscripts (PubMed)](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)

## Data Availability

Unedited Coriell LCL lines are available at https://www.coriell.org/. CRISPR-edited Coriell LCL lines are available upon request and establishment of an MTA with Coriell and NIH/CARD abiding by the Coriell NINDS Human Genetics Repository Material Transfer Agreement For Biospecimens. All generated LCL Coriell ONT DNAseq, CAGEseq and RNAseq data (ILM and ONT) is available at https://www.amp-pd.org/ via GP2 tier 2 access. AMP-PD ILM blood based RNAseq is available at https://www.amp-pd.org/ after signing the data use agreement. 1000 Genomes project data is publicly available at https://www.internationalgenome.org/. All generated brain tissue bulk RNAseq data (ONT) is currently being submitted to the NIMH data sharing platform at https://nda.nih.gov/. Summary statistics for cis-eQTLs and a catalog of ancestry-specific eQTLs from Kachuri et al.12 were obtained from https://doi.org/10.5281/zenodo.7735723. Frontal cortex data publication is in progress and will be available here https://nda.nih.gov/edit_collection.html?id=3151. All of the sequencing data (DOI 10.5281/zenodo.10962119; release 7) used in the preparation of this article were obtained from the Global Parkinson’s Genetics Program (GP2).

## Repository Orientation
- The `analyses/` directory includes all analyses discussed in the manuscript
- The `figures/` directory includes all figures and supplemental figures referenced in the manuscript *(pending publication)*
- The `tables/` directory includes all tables and supplemental tables referenced in the mansucript *(pending publication)*

---

### Analyses

Languages: Bash, R, Python

Each analysis directory has:
1. A notebook that calls the scripts
2. A collection of scripts called in the notebook

An overview of the analysis directories in order:

| Analysis Directory                   | Description                                                                                                                                                                                          |
|----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| OO_ONT_RNA                 | 1) Process raw ONT RNAseq data from basecalling to mapping. 2) Calculate coverage across regions and generate plots. 3) Get TPMs for GBA1 transcripts.                                           |
| 01_Illumina_RNA            | 1) Map Illumina RNA data (fastq). 2) Calculate coverage across regions and generate plots.                                                                                                       |
| 02_ONT_DNA                 | 1) Process raw ONT RNAseq data from base calling to mapping.                                                                                                                                     |
| 03_Branchpoint_predictions | 1) Generate branch point probabilities and generate plot.                                                                                                                                              |
| 04_HBCC_10X                | 1) Process raw ONT single-nuclei data from basecalling to mapping. 2) Split out data by cell types. 3) Calculate coverage across regions and generate plots. 4) Plot barcode sequence diversity. |

---
### Figures and Supplemental Figures

*(pending publication)*

---
### Tables and Supplemental Tables 

*(pending publication)*

---

| Software                            | Version(s)    | Resource URL                                                                                                                                                                                  | RRID            | Notes                                                                                                                           |
|:-----------------------------------:|:-------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:---------------:|:-------------------------------------------------------------------------------------------------------------------------------:|
| bcl2fastq                           | 2.20.0.422    | https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html                                                                                                | RRID:SCR_015058 | Convert Illumina binary files to FASTQ                                                                                          |
| Biorender                           | NA            | biorender.com                                                                                                                                                                                 | RRID:SCR_018361 | Generate Figure 6                                                                                                               |
| Branchpointer                       | 4.3.0         | https://www.bioconductor.org/packages/release/bioc/html/branchpointer.html                                                                                                                    | NA              | Branchpoint prediction algorithm                                                                                                |
| BWA                                 | 0.5.9         | http://bio-bwa.sourceforge.net/                                                                                                                                                               | RRID:SCR_010910 | Mapping  of CAGE-seq data                                                                                                       |
| Clair3                              | 1.0.3         | https://github.com/HKU-BAL/Clair3                                                                                                                                                             | NA              | Small variant calling for long-read DNA sequencing                                                                              |
| DRAGEN                              | 3.7.8         | https://developer.illumina.com/dragen                                                                                                                                                         | NA              | Variant caller for Illumina data                                                                                                |
| ggplot2                             | 3.5.1         | https://cran.r-project.org/web/packages/ggplot2/index.html                                                                                                                                    | RRID:SCR_014601 | Generation of R figures                                                                                                         |
| gnomAD                              | 4.0           | http://gnomad.broadinstitute.org/                                                                                                                                                             | RRID:SCR_014964 | Exploration of reported GBA1 variants                                                                                           |
| Guppy                               | 6.1.2         | https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview                                                   | RRID:SCR_023196 | Basecalling Oxford Nanopore R9.4.1 files                                                                                        |
| Integrative Genome Browser          | 2.16.0        | https://igv.org/                                                                                                                                                                              | RRID:SCR_011793 | Visualization of mapped sequenced files and transcripts                                                                         |
| Mascot                              | NA            | http://www.matrixscience.com/server.html                                                                                                                                                      | RRID:SCR_014322 | Protein search engine                                                                                                           |
| Minimap2                            | 2.24 and 2.26 | https://github.com/lh3/minimap2                                                                                                                                                               | RRID:SCR_018550 | Mapping of long-read sequencing files                                                                                           |
| Minknow                             | 22.10.7       | https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revao_14dec2018                                                                           | RRID:SCR_023196 | Oxford Nanopore sequencing software                                                                                             |
| Plink                               | 2.0           | http://www.nitrc.org/projects/plink                                                                                                                                                           | RRID:SCR_001757 | Extracting rs3115534 genotypes                                                                                                  |
| Proteome Discoverer                 | 2.4           | https://www.thermofisher.com/us/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html | RRID:SCR_014477 | Protein database search                                                                                                         |
| Pychopper                           | 2.7.1         | https://github.com/epi2me-labs/pychopper                                                                                                                                                      | RRID:SCR_018966 | Identify, orient, and trim Oxford Nanopore cDNA full length reads                                                               |
| PycoQC                              | 2.5.2         | https://a-slide.github.io/pycoQC/                                                                                                                                                             | RRID:SCR_024185 | Generate quality reports for Oxford Nanopore sequencing runs                                                                    |
| Python Programming Language         | 3.1           | http://www.python.org/                                                                                                                                                                        | RRID:SCR_008394 | pandas;seqio; Used for sorting through fastq files and data manipulation                                                        |
| R Project for Statistical Computing | 4.3           | http://www.r-project.org/                                                                                                                                                                     | RRID:SCR_001905 | branchpointer;dplyr;ggplot2;Used for generating coverage plots, running regression, and calculating branch point probabilities. |
| Real Time Analysis                  | 3.0           | http://support.illumina.com/sequencing/sequencing_software/real-time_analysis_rta.html                                                                                                        | RRID:SCR_014332 | Basecall raw images generated by Illumina                                                                                       |
| RegulomeDB                          | 2.2           | http://www.regulomedb.org/                                                                                                                                                                    | RRID:SCR_017905 | Assess elements in non-coding regions                                                                                           |
| Samtools                            | 1.17          | http://www.htslib.org/                                                                                                                                                                        | RRID:SCR_002105 | Calculate sequencing depth, subset bams, convert bams to fastq                                                                  |
| Seqkit                              | 2.2.0         | https://bioinf.shenwei.me/seqkit/                                                                                                                                                             | RRID:SCR_018926 | Generate sequencing statistics from Oxford Nanopore raw data                                                                    |
| Sniffles                            | 2.2           | https://github.com/fritzsedlazeck/Sniffles                                                                                                                                                    | RRID:SCR_017619 | Call structural variants from long-read data                                                                                    |
| SpliceAI                            | 1.3           | https://pypi.org/project/spliceai/                                                                                                                                                            | RRID:SCR_024532 | Predict splice junctions in mRNA                                                                                                |
| STAR                                | 2.7.10        | https://github.com/alexdobin/STAR                                                                                                                                                             | RRID:SCR_004463 | Mapping of short-read sequencing files                                                                                          |
| Stringtie2                          | 2.2.1         | https://ccb.jhu.edu/software/stringtie/                                                                                                                                                       | RRID:SCR_016323 | Transcript calling and quantification                                                                                           |
| UCSC Genome Browser                 |               |                                                                                                                                                                                               |                 |                                                                                                                                 |
