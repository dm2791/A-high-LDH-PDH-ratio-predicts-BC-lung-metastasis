# Divergent-use-of-metabolic-fluxes-in-breast-cancer-metastasis
MATLAB scripts to analyze the data and make the figures for the paper "Divergent use of metabolic fluxes in breast cancer metastasis".

This repository includes some of the data analyzed. 

figure1_metabolomics_and_mithril/breastnorm.xlsx: Metabolomics data produced by Metabolon.
figure1_metabolomics_and_mithril/geneset.txt: the genes in the GLYCOLYSIS gene set
figure1_metabolomics_and_mithril/mithriloutputbrainh.xls: output file from mithril analysis
figure1_metabolomics_and_mithril/mithriloutputbrm2lm2.txt: output file from mithril analysis
figure1_metabolomics_and_mithril/output_mithril_lung.txt: output file from mithril analysis
figure2_fluxes_measured/231_mito.xlsx: results from seahorse mitochondria assay
figure2_fluxes_measured/Deepti_231_test_GlycolyticRate_050219_norm.xlsx: results from seahorse glycolysis assay
figure2_fluxes_measured/YSI_data.xlsx: results from YSI assay
figure3_fba_enzyme_activities/Flexibility-XF Mito Fuel Flex Test_231test2_071019.xlsx: results from seahorse
figure3_fba_enzyme_activities/G3PDHREPEAT0MIN.xlsx: G3PDH enzyme activities 
figure3_fba_enzyme_activities/PK0min.xlsx: PK enzyme activities 
figure3_fba_enzyme_activities/PFK5min.xlsx: PFK enzyme activities 
figure3_fba_enzyme_activities/hexokinase50min.xlsx: hexokinase enzyme activities 
figure4_LDHPDH_ratio/data_clinical_patient.txt: Clinical metadata from Metastatic Breast Cancer Project

The following large files must be downloaded by the user:

1) For codes involving mRNA analysis: affymetrix chip annotation
Download from: http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu133
HG-U133A Annotations, CSV format, Release 36 (19 MB, 7/12/16)
Save as .txt file

2) For Metastatic Breast Cancer Project Analysis: patient data files
https://www.cbioportal.org/study/summary?id=brca_mbcproject_wagle_2017
(MBCproject cBioPortal data version March 2020)


