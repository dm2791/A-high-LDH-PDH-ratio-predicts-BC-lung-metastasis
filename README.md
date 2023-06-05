# A-high-LDH-PDH-ratio-predicts-BC-lung-metastasis
MATLAB scripts to analyze the data and make the figures for the paper "The ratio of key metabolic transcripts is a predictive biomarker of breast cancer metastasis to the lung", running title "A high LDH/PDH ratio predicts BC lung metastasis"

Associated with this repository is the following metabolomics data release:
https://zenodo.org/record/7996343

This repository includes some of the data analyzed: 

AACRprojectGENIE/mbc_genie_2020_clinical_data.csv: data from AACR project GENIE
AACRprojectGENIE/Mutated_Genes.txt: mutation data from AACR project GENIE

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

enzyme_activities/g3pdhdatamat.csv: G3PDH enzyme activities
enzyme_activities/pkdatamat.csv: PK enzyme activities
enzyme_activities/pfkdatamat.csv: PFK enzyme activities
enzyme_activities/hexokinasedatamat.csv: hexokinase enzyme activities

MBCP_analysis_and_LDHPDH_ratio/ldhbrain.csv: transcriptomic output
MBCP_analysis_and_LDHPDH_ratio/ldhlung.csv: transcriptomic output
MBCP_analysis_and_LDHPDH_ratio/pdhbrain.csv: transcriptomic output
MBCP_analysis_and_LDHPDH_ratio/pdhlung.csv: transcriptomic output

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
