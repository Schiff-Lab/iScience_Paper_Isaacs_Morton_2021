# STEP 1.0: INSTALLATION

## STEP 1.1: Install Bioconductor and GSVA/Limma tools
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase", "GSVA", "Mfuzz", "NMF", "gage", "limma"))

## STEP 1.2: Install ProteoQ packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("qzhang503/proteoQ")


# STEP 2.0: EXPERIMENTAL SETUP, DATA PREPROCESSING AND QA

## STEP 2.1: Set up the experiments

### STEP 2.1.1 Load proteoQ library
library(proteoQ)

### STEP 2.1.2 Specifiy the working directory, which should contain the following:
####  1. PSM table (csv format) which was obtained from the Mascot search engine [It is expected the filename will be in Mascot format i.e. starting with the letter F and followed by 6 digits e.g. F003499.csv]
####  2. Metadata of samples (excel or csv) which must contain a) multiplex experiment numbers, b) TMT channels, c) LC/MS injection indices, d) sample IDs, e) Corresponding RAW data filenames, f) any other metadata such as clinical information etc.
####  3. Fractionation data (excel or csv) --> This file is optional and only required if the samples were fractionated offline prior to the LC/MS run. Type ?setup_expts for more details
dat_dir <- "C:\\Users\\aisaac01\\Desktop\\CONSHA_proteoQ"
system.file("extdata","F003499.csv", package = "proteoQ")
system.file("extdata","expt_smry.xlsx", package = "proteoQ")
system.file("extdata","frac_smry.xlsx", package = "proteoQ")

### STEP 2.1.3 Load the experiment: Should expect c2_msig_hs, dbs, fraction_scheme, label_scheme, label_scheme_full to load in the environemt if using R Studio
setup_expts()

## STEP 2.2 Summarize PSMs to peptides
### STEP 2.2.1 Generate PSM reports: Should expect the PSM data to be normalized
normPSM(
  rptr_intco = 1000, 
  rm_craps = FALSE, 
  rm_krts = FALSE, 
  rm_outliers = FALSE, 
  plot_violins = TRUE
)

## STEP 2.2.2 Generate peptide reports: Should expect the peptide data to be normalized
normPep(
  id = pep_seq_mod, 
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(10, 95), 
  range_int = c(5, 95), 
  n_comp = 2, 
  seed = 321, 
  annot_kinases = TRUE, 
  maxit = 200, 
  epsilon = 1e-05
)

## STEP 2.2.3 Generate 2 plots, one with and the other without scaling of Log2FC to decide if scaling of peptide or PSM data is necessary for the dataset

### STEP 2.2.3.1 Plot without scaling of Log2FC
pepHist(
  scale_log2r = FALSE,
  show_curves = TRUE,
  show_vline = TRUE,
  ncol = 5
)

### STEP 2.2.3.2 Plot with scaling of Log2FC
pepHist(
  scale_log2r = TRUE,
  show_curves = TRUE,
  show_vline = TRUE,
  ncol = 5
)

## STEP 2.3 Summarize peptides to proteins

### STEP 2.3.1 Generate protein reports
normPrn(
  id = gene,
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 90),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 246,
  fasta = "C:\\Results\\DB\\Refseq\\RefSeq_HM_Frozen_20130727.fasta",
  maxit = 200,
  epsilon = 1e-05
)

## STEP 2.3.2 Generate 2 plots, one with and the other without scaling of Log2FC to decide if scaling of protein data is necessary for the dataset

### STEP 2.3.2.1 Plot without scaling of Log2FC
prnHist(
  scale_log2r = FALSE,
  show_curves = TRUE,
  show_vline = TRUE,
  ncol = 5
)

### STEP 2.3.2.1 Plot with scaling of Log2FC
prnHist(
  scale_log2r = TTRUE,
  show_curves = TRUE,
  show_vline = TRUE,
  ncol = 5
)

## STEP 2.4: Correlate intensity with Log2FC for both peptides and proteins

### STEP 2.4.1: Correlation plots of peptide data
pepCorr(
  use_log10 = TRUE,
  scale_log2r = TRUE,
  min_int = 3.5,
  max_int = 6.5,
  min_log2r = -2,
  max_log2r = 2,
  width = 24,
  height = 24
)

### STEP 2.4.2: Correlation plots of protein data
prnCorr(
  use_log10 = TRUE,
  scale_log2r = TRUE,
  min_int = 3.5,
  max_int = 6.5,
  min_log2r = -2,
  max_log2r = 2,
  width = 24,
  height = 24
)


# STEP 3.0: PEPTIDE ANALYSIS

## STEP 3.1: MDS plots of peptides
pepMDS(
  scale_log2r = TRUE, 
  adjEucDist = FALSE, 
  show_ids = FALSE
)

## STEP 3.2: PCA plots of peptides
pepPCA(
  scale_log2r = TRUE, 
  show_ids = FALSE
)

## STEP 3.3: Eucledian distance of peptides
pepEucDist(
  scale_log2r = TRUE, 
  adjEucDist = FALSE, 
  show_ids = FALSE, 
  annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
  
  display_numbers = FALSE, 
  number_color = "grey30", 
  number_format = "%.2f",
  
  clustering_distance_rows = "euclidean", 
  clustering_distance_cols = "euclidean", 
  
  fontsize = 16, 
  fontsize_row = 20, 
  fontsize_col = 20, 
  fontsize_number = 16, 
  
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "grey60", 
  cellwidth = 20, 
  cellheight = 20, 
  width = 40, 
  height = 40
)

## STEP 3.4: Peptide heatmap
pepHM(
  scale_log2r = TRUE, 
  
  xmin = -2, 
  xmax = 2, 
  x_margin = 0.2, 
  
  annot_cols = c("Peptide_Yield", "TMT_Set", "Group"), 
  cluster_rows = TRUE, 
  cutree_rows = 6, 
  show_rownames = FALSE, 
  show_colnames = TRUE, 
  fontsize_row = 3, 
  cellwidth = 14, 
  width = 24, 
  height = 12
)


# STEP 4.0: PROTEIN ANALYSIS

## STEP 4.1: MDS plots of proteins
prnMDS(
  scale_log2r = TRUE, 
  adjEucDist = FALSE, 
  show_ids = FALSE
)

## STEP 4.2: PCA plots of proteins
prnPCA(
  scale_log2r = TRUE, 
  show_ids = FALSE
)

## STEP 4.3: Eucledian distance of proteins
prnEucDist(
  scale_log2r = TRUE, 
  adjEucDist = FALSE, 
  show_ids = FALSE, 
  annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
  
  display_numbers = FALSE, 
  number_color = "grey30", 
  number_format = "%.2f",
  
  clustering_distance_rows = "euclidean", 
  clustering_distance_cols = "euclidean", 
  
  fontsize = 16, 
  fontsize_row = 20, 
  fontsize_col = 20, 
  fontsize_number = 16, 
  
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "grey60", 
  cellwidth = 20, 
  cellheight = 20, 
  width = 40, 
  height = 40
)

## STEP 4.4 Protein heatmap
prnHM(
  scale_log2r = TRUE, 
  xmin = -2, 
  xmax = 2, 
  x_margin = 0.2, 
  annot_cols = c("Peptide_Yield", "TMT_Set", "Group"), 
  cluster_rows = TRUE, 
  cutree_rows = 6, 
  show_rownames = FALSE, 
  show_colnames = TRUE, 
  fontsize_row = 3, 
  cellwidth = 14, 
  width = 24, 
  height = 12
)

# STEP 5.0: IMPUTATION

## STEP 5.1: Peptide data imputation
pepImp(m = 5, maxit = 5)

## STEP 5.1: Protein data imputation
prnImp(m = 5, maxit = 5)

# STEP 6.0: CALCULATE P VALUES BETWEEN GROUPS (PROTEIN LEVEL)
prnSig(
  scale_log2r = TRUE, 
  formula_1 = ~ Term["Paeni_PIH-NPIH", "NonPaeni_PIH-NPIH", "Paeni_PIH-NonPaeni_PIH"]
)

# STEP 7.0: CREATE VOLCANO PLOTS (PROTEIN LEVEL)
prnVol(scale_log2r = TRUE)

# STEP 8: PLOT TREND CLUSTERS
prnTrend(n_clust = 6, scale_log2r = TRUE)

# STEP 9: PERFORM NONNEGATIVE MATRIX FACTORIZATION (NMF)
if (!"NMF" %in% installed.packages()[, "Package"]) install.packages("NMF")
library(NMF)
prnNMF(r = 6, xmin = -2, xmax = 2, x_margin = 0.2, 
       annot_cols = c("Peptide_Yield", "TMT_Set", "Group"), 
       scale_log2r = TRUE)

# STEP 10: RUN AND PLOT GENE SET VARIATION ANALYSIS (USING GO, KEGG AND C2MSIG)
prnGSVA(gset_nm = c("go_sets", "kegg_sets", "c2_msig"), scale_log2r = TRUE)
gsvaMap(scale_log2r = TRUE, pval_cutoff = 1E-8, show_sig = "pVal")

# STEP 11: RUN AND PLOT GENE SET ENRICHMENT (GAGE) ANALYSIS (USING GO, KEGG AND C2MSIG)
prnGAGE(gset_nm = c("go_sets", "kegg_sets", "c2_msig"), scale_log2r = TRUE)
gageMap(scale_log2r = TRUE)
