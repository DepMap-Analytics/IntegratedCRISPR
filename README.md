# IntegratedCRISPR
Code for integration and evaluation of Broad and Sanger Institutes CRISPR-Cas9 screens 


The code should be executed in the following order:
GenerateSourceData.R
Figure2.R
Figure3.R
UseCaseEssNonEss.R
CalculateNormLRT.R
UseCaseLineageSubtypes.R
UseCaseBiomarkers.R
CalculateBinaryDepletions.R


To run the code, source data files should be downloaded from the Figshare repository: 
The location of the Figshare directory should be added to the R scripts as indicated (dir.Input variable).

The packages required to run the R scripts are:
here
preprocessCore
stringr
RColorBrewer
colorspace
scales
psych
sva
ggfortify
VennDiagram
splines
e1071
tsne
caTools
ggplot2
gridExtra
magrittr
tidyverse
dplyr
MASS
sn
