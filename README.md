# IntegratedCRISPR
Code for integration and evaluation of Broad and Sanger Institutes CRISPR-Cas9 screens. These scripts are to reproduce the analysis in the paper. For the integrations functions alone there is an R package available here: https://github.com/DepMap-Analytics/IntCRISPR


The code should be executed in the following order:
* GenerateSourceData.R
* Figure2.R
* Figure3.R
* UseCaseEssNonEss.R
* CalculateNormLRT.R
* UseCaseLineageSubtypes.R
* UseCaseBiomarkers.R
* CalculateBinaryDepletions.R
* CalcCommonEssentials.R
* IntegratedPerformance.R
* DownSamplingBatchCorrection.R
* DownSamplingFigures.R

To run the code, source data files should be downloaded from the Figshare repository: 
https://figshare.com/projects/Integrated_CRISPR/78252
The location of the Figshare directory should be added to the R scripts as indicated (dir.Input variable). Note the DownSamplingBatchCorrection script is set to run in parallel over 24 cores. This script is computationally expensive.

All python scripts are contained in IntegrationAnalysis.ipynb. This file requires jupyter notebook to open. All the data inputs required to run the Jupiter notebook are in the Figshare (https://figshare.com/projects/Integrated_CRISPR/78252) subdirectory NotebookInput.

The packages required to run the R scripts are:
* here
* preprocessCore
* stringr
* RColorBrewer
* colorspace
* scales
* psych
* sva
* ggfortify
* VennDiagram
* splines
* e1071
* tsne
* caTools
* ggplot2
* gridExtra
* magrittr
* tidyverse
* dplyr
* MASS
* sn
* beeswarm
