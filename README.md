# NRF2-drug-repositioning
Computational pipeline for identifying NRF2 inhibitors using transcriptomic signature reversal

This repository contains the scripts used for computational drug repositioning targeting NRF2 in hepatocellular carcinoma.

## Contents
- NRF2 consensus signature construction (TCGA + LINCS)
- Connectivity-based drug repositioning analysis
- CKAP4 validation analysis

## Requirements
- R (version 4.5.2)
- Packages: dplyr, data.table, clusterProfiler, cmapR

## Generalisability
Although developed for NRF2, this pipeline can be applied to other gene expression signatures, allowing identification of compounds that reverse alternative transcriptional programmes.

## Notes
Large datasets (e.g. LINCS L1000) are not included and must be downloaded separately.
