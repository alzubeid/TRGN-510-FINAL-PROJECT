# TRGN-510-FINAL-PROJECT

# Title

Genomic and transcriptomic analysis of BRCA-associated triple negative breast cancer.

# Author:

Batul Al-zubeidy

# Overview:

While triple negative breast cancer accounts for about 15% of total invasive breast cancer cases diagnosed in the US every year, it represents approximately 30% of BRCA-associated breast tumors.  This aggressive subtype is characterized by lack of expression of the estrogen, and progesterone hormone receptors (both at <1% expression), as well as decreased expression of the membrane surface protein Human Epidermal Growth Factor Receptor 2 (HER2) on epithelial mammary cells.    Triple negative tumors carry worse prognosis compared with hormone-receptor positive cancers with 5-year-survival rates ranging from 12%-91% for stage I-IV, respectively (American Cancer Society).  
https://www.cancer.org/cancer/breast-cancer/about/types-of-breast-cancer/triple-negative.html
	These tumor islands express significant heterogeneity reflecting various genetic variations that likely explain their variable clinical response/outcome.  This project will utilize TCGA high throughput genomic data to investigate BRCA female patients (n=10) diagnosed with triple negative breast cancer between the ages of 26 and 90.  Similarily, I will analyze BRCA female patients diagnosed with ER+ PR+ HER2-non-amplified (n=10) diagnosed between the ages of 26 and 90, too.  I will be controlling for race and ethnicity.  
I will be analyzing gene expression through RNA sequencing to identify the mutational burdens as possible predictors of clinical outcome (such as stage, clinical response, and disease free survival).  

# Data:

Use cBioPortal combined studies: 
-Breast invasive carcinoma (TCGA, Cell 2015)
-Breast invasive carcinoma (TCGA, Firehose legacy)
-Breast invasive carcinoma (TCGA, Nature 2012)
-Breast invasive carcinoma (TCGA, PanCancer Atlas)

-CSV file with the included triple-negative breast cancer samples is attached as "combined_study_clinical_data-3.csv"
-CSV file with the included ER+ PR+ HER2-non-amplified breast cancer samples is attached as "combined_study_clinical_data.csv"

Clinical data obtained from (cBioPortal):
https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pub2015%2Cbrca_tcga%2Cbrca_tcga_pub%2Cbrca_tcga_pan_can_atlas_2018

Genomic/transcriptomics from: 
https://portal.gdc.cancer.gov

I will utilize the package DeSEQ2 (http://bioconductor.org/packages/release/bioc/vignettes/DESeq2.html)

# Milestone #1:
-Due to the large clinical dataset involving the two subsets of patients, a smaller subset was chosen to run a test relative to the vignette.  
-Differential expression analysis was run on the small selected sample
	-a subset of 10 patients with ER+/PR+/HER2- and 10 patients with triple negative breast cancer (the 2 attached file are TNB_HR_Joined_Counts test.csv, and metadata test.csv).  

# Milestone #2:
-Differential expression analysis was run on the small selected sample
	-a subset of 10 patients with ER+/PR+/HER2- and 10 patients with triple negative breast cancer: 
	showing 1,304 genes  with stastically significant difference.  
	-PCA plot show that the TN breast cancers express distinct genes from hormone-receptor-positive tumors.  
	
	
# Conclusion:
1. Triple-negative breast cancer in females appeared to show significant differences in gene expression relative to hormone-receptor-positive breast cancers. 


2. The principal component separation plot showed that the two cohorts significantly different. 
3. These differences may overestimate the true genetic alterations between the histological tissue.  This overestimation may be due to racial, age-related, or hormonal status.
4. Further study including larger sample sizes should be performed, where race, hormonal status, and stage can be compared. 


# Known Issues:
-Small sample sizes were used to run the test for differential expression, which likely obsecures the true gene expression.  
-The data obtained from TCGA cases were cross examined with portal.gdc.cancer.gov dataset to obtain the clinical pathology, including estrogn/progestrone/HER2 biomarkers.  
-The RNA-sequence data showed significantly different PCA analyses between triple negative and estorgen-receptor-positive breast cancer; much greater differences than anticipated.  

Requirements:
This project requires RStudio.
The data was downloaded from cBioPortal (https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pub2015%2Cbrca_tcga%2Cbrca_tcga_pub%2Cbrca_tcga_pan_can_atlas_2018) and Genomic/transcriptomics from (https://portal.gdc.cancer.gov). 
There are various installations required (see Final project for the entire list), which were downloaded from http://bioconductor.org/packages
