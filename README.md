# Integrative Multi-Omic Analysis Using Coupled Matrix Factorization (CMF)

This repository contains the code and resources for performing integrative multi-omic analysis using coupled matrix factorization (CMF) to study the molecular mechanisms of arsenic toxicity. 

## Overview
Arsenic exposure disrupts the epigenome, transcriptome, and metabolome, leading to severe health implications. This project leverages the **PARAFAC2 AO-ADMM model** for CMF to integrate and analyze:
- **RRBS**: Reduced Representation Bisulfite Sequencing data.
- **RNA-seq**: Transcriptomic data.
- **Metabolomics**: Metabolic profiling data.

By jointly factorizing these datasets, we uncover regulatory mechanisms driving molecular changes in mouse ESCs and EpiLCs treated with arsenic.

## Problem Statement
The study aims to:
1. Identify methylation dysregulation trends that influence transcriptional and metabolic changes.
2. Reveal distinct regulatory networks linked to arsenic-induced toxicity and its pathophysiological effects, including cancer and developmental disorders.

## Methods
We implemented CMF using the **PARAFAC2 AO-ADMM algorithm**, which allows:
- Joint factorization of datasets sharing common features but differing in row dimensions.
- Regularization for robust decomposition and interpretability of factor matrices.

The optimization problem is solved with l1 regularization to ensure stability and uniqueness in the factorized components.
