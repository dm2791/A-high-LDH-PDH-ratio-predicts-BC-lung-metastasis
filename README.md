# Divergent-use-of-metabolic-fluxes-in-breast-cancer-metastasis
MATLAB scripts to analyze the data and make the figures for the paper "Divergent use of metabolic fluxes in breast cancer metastasis".

This repository includes some of the data analyzed: 

figure1_metabolomics_and_mithril/breast.csv: Metabolomics data produced by Metabolon.
figure1_metabolomics_and_mithril/hallmarkgeneset.csv: the genes in the HALLMARK GLYCOLYSIS gene set
figure1_metabolomics_and_mithril/brainoutput_mithril.csv: output file from mithril analysis
figure1_metabolomics_and_mithril/lungoutput_mithril.csv: output file from mithril analysis
figure1_metabolomics_and_mithril/brm2lm2output_mithril.csv: output file from mithril analysis
figure1_metabolomics_and_mithril/brm2gsea.csv: transcriptomic data, brm2 vs parental
figure1_metabolomics_and_mithril/lm2gsea.csv: transcriptomic data, lm2 vs parental
figure2_fluxes_measured/mitodata_O2.csv: results from seahorse mitochondria assay
figure2_fluxes_measured/mitodata_ECAR.csv: results from seahorse mitochondria assay
figure2_fluxes_measured/glycolyticratedata_O2.csv: results from seahorse glycolysis assay
figure2_fluxes_measured/glycolyticratedata_PER.csv: results from seahorse glycolysis assay
figure2_fluxes_measured/YSIdatamat.csv: results from YSI assay
figure3_fba_enzyme_activities/capacitymat.csv: results from seahorse
figure3_fba_enzyme_activities/dependencymat.csv: results from seahorse
figure3_fba_enzyme_activities/g3pdhdatamat.csv: G3PDH enzyme activities 
figure3_fba_enzyme_activities/pkdatamat.csv: PK enzyme activities 
figure3_fba_enzyme_activities/pfkdatamat.csv: PFK enzyme activities 
figure3_fba_enzyme_activities/hexokinasedatamat.csv: hexokinase enzyme activities 
figure4_LDHPDH_ratio/ldhbrain.csv: transcriptomic output
figure4_LDHPDH_ratio/ldhlung.csv: transcriptomic output
figure4_LDHPDH_ratio/pdhbrain.csv: transcriptomic output
figure4_LDHPDH_ratio/pdhlung.csv: transcriptomic output

The following large files must be downloaded by the user:

1) For codes involving mRNA analysis: affymetrix chip annotation
Download from: http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu133
HG-U133A Annotations, CSV format, Release 36 (19 MB, 7/12/16)
Save as .txt file

2) For Metastatic Breast Cancer Project Analysis: patient data files
https://www.cbioportal.org/study/summary?id=brca_mbcproject_wagle_2017
(MBCproject cBioPortal data version March 2020)

3) For Metastatic Breast Cancer Project Analysis: UMAP code
https://umap-learn.readthedocs.io/en/latest/


Computational analysis of metabolic fluxes was implmented in the Python notebook:
figure3_fba_enzyme_activities/fba/main.ipynb
The directory figure3_fba_enzyme_activities/fba/ is a standalone folder that includes all necessary input files for
the Python notebook to run. These input files are: (1) metabolic model (mammalian_cell_model_final.xlsx); (2) YSI data (YSI.xlsx);
(3) Seahorse data (Seahorse_atp.xlsx). Installation of the following python packages is required: pandas, numpy, seaborn, cobra, sklearn.
