# Git_Differentiation_Metabolism
# README : Bioinformatic analysis of the project

********************************************************************************
********************************************************************************

## I. Folder tree

Follows the working structure :   
- bin folder : contains the scripts (R, Snakemake...)  
- data folder : contains raw data  
- report folder : contains html report from each script of the bin folder    
- exp folder : contains script's outputs (plots, tables...)    
- source codes to launch all scripts in the correct order    
- one or more README files  
  
Data and bin folders contain subfolders from each data type : scATACseq, scRNAseq or bulk ATACseq.  
Raw files from data folder were renamed for clarity (name of the condition in each file name).  
Exp folder contains subfolders named after their associated script. In each subfolders are sub-subfolders named according to their date of creation.  
Html files are dated with the date on which they were created.  

## II. Raw data description

This project contains multi-omics data according to different metabolic stress conditions.   

### scATAC-seq 10X

Experiment done on the 25th of October 2021  
**4 conditions** : CTRL, DON, 2DG, AOA    
**Timing** : 24h    
.fastq files were treated by Imagine's bioinformatic plateform.   
Our analysis here begins with the output files from Cell Ranger and is carried out with R.  
  
Initially, data were collected from Cell Ranger's output folders named full_cellranger_outs/ATAC_DrugName. For simplicity, we copy-paste in the git's data folder only useful files for further analysis. In addition, drug's name were added to the beggining of the file names. Finally, barcodes.tsv, matrix.mtx and peaks.bed files that were stored in a filtered_peak_bc_matrix folder were taken out and the subfolder's name was added to the name of each individual file.   


### scRNA-seq 10X with CITE-seq (CD34, CD133)

Experiments done on the 27th of May 27 2021 and 04th of July 04 2022  
**10 conditions** : CTRL x2, DON, 2DG, AOA, CTRLaK, DONaK, 2DGaK, AOAaK, VPA      
**Timing** : 96h      
.fastq files were treated by Cochin's bioinformatic plateform GENOM'IC.   
Our analysis here begins with the output files from Cell Ranger and is carried out with R.  
  

### bulk ATAC-seq with FAST-ATAC method

Experiments done on the 23th of July 23 2019 and the 29th of October 29 2019.  
**7 conditions** : Xvivo, CTRL, DON, 2DG, AOA, CTRLaK, VPA    
**Timings**: 00h, 03h, 06h, 12h, 24h      
.fastq files were treated by CEA's bioinformatic plateform.     
Our analysis here begins with .bam files and is carried out with Snakemake and R.  
  

## III. Analysis description

Refer to README_workflow.pptx to understand script's order.  
Read src files and use them to launch all the analysis.  
