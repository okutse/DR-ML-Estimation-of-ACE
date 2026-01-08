# Title: Doubly-robust machine learning estimation of the average causal effect of (selected) DEGs on breast cancer stage and survival outcomes
# Author: Amos Okutse
# Date: Jan 2025


# Comparative Methods:
#(1) DR Penalized Logistic Regression
#(2) DR Random Forest
#(3) DR BART
#(4) DR Logistic regression on categorized DEGs (low, mid, high expression levels based on quantiles)

## Load required packages
library(data.table)    # fast I/O and merging
library(dplyr)         # data manipulation
library(stringr)       # string helpers
library(survival)      # survival objects
library(glmnet)        # high-dimensional regression
library(grf)          # generalized random forests
library(BART)         # Bayesian Additive Regression Trees
library(splines)      # natural splines
library(ggplot2)      # plotting
library(ggrepel)      # labeling points cleanly
library(ggpubr)       # arranging plots
library(cowplot)      # arranging plots
library(parallel)     # parallel computing
library(doParallel)   # parallel backend for foreach
library(foreach)      # parallel loops
library(boot)        # for inverse logit function

# install BiocManager if you don't have it (uncomment and run only run once)
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# install the limma package using BiocManager
# BiocManager::install("limma")
library(limma)       # for DEG analysis
# library(edgeR)       # for DEG analysis on NB counts
library(matrixStats) # for simple matrix operations
library(VennDiagram)  # for Venn diagrams
library(glmnet)       # high-dimensional regression


# source helper files
source("scripts/helper.R")

## Load the TCGA BRCA data files from data directory in DataFest root directory
cp  <- data.table::fread("data/Data_clinical_patient.txt")[5:.N] # drop 1-4th rows
cs   <- data.table::fread("data/Data_clinical_sample.txt")[5:.N] # drop 1-4th rows
mp     <- fread("data/Data_methylation_promoters_rrbs.txt") # maybe augment with PCA exp (90%)
mexp     <- data.table::fread("data/Data_mrna_illumina_microarray.txt")
mt      <- data.table::fread("data/data_mutations.txt")                         


################################################################################
# Clinical patient and sample data ----
################################################################################

# check any missing values in cp columns
# Find all "" and replace with NA
cs <- cs %>%
  dplyr::mutate_all(~dplyr::na_if(., ""))

## Rename the cs 
cs <- cs %>% dplyr::mutate(`HER2 Status` = dplyr::na_if(`HER2 Status`, "")) %>% 
  dplyr::rename(patient_id = `#Patient Identifier`,
              sample = `Sample Identifier`,
              cancer_type = `Cancer Type`,
              cancer_type_detailed = `Cancer Type Detailed`,
              er_status = `ER Status`,
              her2_status = `HER2 Status`,
              neoplasm_histologic_grade = `Neoplasm Histologic Grade`,
              oncotree_code = `Oncotree Code`,
              pr_status = `PR Status`,
              sample_type = `Sample Type`,
              tumor_size = `Tumor Size`,
              tumor_stg = `Tumor Stage`,
              tmb_nonsynonymous = `TMB (nonsynonymous)`) %>% 
  # add binarized tumor stage from levels levels(as.factor(cs$`Tumor Stage`)) "0" "1" "2" "3" "4"
  dplyr::mutate(tumor_stage = as.factor(case_when(
    str_detect(tumor_stg, "0|1|2")  ~ "early",
    str_detect(tumor_stg, "3|4") ~ "late",
    TRUE ~ NA_character_
  ))) %>% 
  # add observability mask for tumor stage (1 = observed, 0 = missing) to address missing outcome data
  dplyr::mutate(tumor_stage_mask = if_else(!is.na(tumor_stage), 1, 0)) %>% 
  # add inverse weights for tumor stage to address class imbalance
  dplyr::mutate(tumor_stage_invwt = case_when(
    tumor_stage == "early" ~ 1 / sum(tumor_stage == "early", na.rm = TRUE),
    tumor_stage == "late" ~ 1 / sum(tumor_stage == "late", na.rm = TRUE),
    TRUE ~ NA_real_
  )) %>% 
  # add ER status mask (1 = observed, 0 = missing) to address missing confounder data
  dplyr::mutate(er_status_mask = if_else(!is.na(er_status), 1, 0)) %>% 
  # convert covariates to factors
  dplyr::mutate(
    #cancer_type = as.factor(cancer_type), # all breast cancer hence not meaningful
    cancer_type_detailed = as.factor(cancer_type_detailed),
    er_status = as.factor(er_status),
    her2_status = as.factor(her2_status), 
    her2_status_mask = if_else(!is.na(her2_status), 1, 0),
    neoplasm_histologic_grade = as.factor(neoplasm_histologic_grade),
    neoplasm_histologic_grade_mask = if_else(!is.na(neoplasm_histologic_grade), 1, 0),
    oncotree_code = as.factor(oncotree_code), # no missing
    pr_status = as.factor(pr_status),
    pr_status_mask = if_else(!is.na(pr_status), 1, 0),
    # sample_type = as.factor(sample_type) # all primary hence not meaningful
    tumor_size = as.numeric(tumor_size),
    tumor_size_mask = if_else(!is.na(tumor_size), 1, 0),
    tumor_stg = as.numeric(tumor_stg)
  ) %>% 
  dplyr::select(-c(cancer_type, sample_type)) # drop non-informative covariates

sapply(cs , function(x) sum(is.na(x)) )

# rename the cp data file
# Find all "" and replace with NA
cp <- cp %>%
  dplyr::mutate_all(~dplyr::na_if(., ""))

cp <- cp %>%
  #dplyr::mutate(`Sex` = dplyr::na_if(`Sex`, "")) %>% all female hence not meaningful
  dplyr::rename(patient_id = `#Patient Identifier`,
              lymph_nodes_examined_positive = `Lymph nodes examined positive`,
              nottingham_prognostic_index = `Nottingham prognostic index`,
              cellularity = `Cellularity`,
              chemotherapy = `Chemotherapy`,
              cohort = `Cohort`,
              er_status_measured_by_ihc = `ER status measured by IHC`,
              her2_status_measured_by_snp6 = `HER2 status measured by SNP6`,
              hormone_therapy = `Hormone Therapy`,
              inferred_menopausal_state = `Inferred Menopausal State`,
              #sex = `Sex`, # not meaningful
              integrative_cluster = `Integrative Cluster`,
              age_at_diagnosis = `Age at Diagnosis`,
              overall_survival_months = `Overall Survival (Months)`,
              overall_survival_status = `Overall Survival Status`,
              pam50_claudin_low_subtype = `Pam50 + Claudin-low subtype`,
              three_gene_classifier_subtype = `3-Gene classifier subtype`,
              patients_vital_status = `Patient's Vital Status`,
              primary_tumor_laterality = `Primary Tumor Laterality`,
              radio_therapy = `Radio Therapy`,
              tumor_other_histologic_subtype = `Tumor Other Histologic Subtype`,
              type_of_breast_surgery = `Type of Breast Surgery`,
              relapse_free_status_months = `Relapse Free Status (Months)`,
              relapse_free_status = `Relapse Free Status`) %>% 
dplyr::mutate(
              # categorical features
              relapse_free_status = as.factor(relapse_free_status),
              overall_survival_status = as.factor(overall_survival_status),
              cellularity = as.factor(cellularity),
              chemotherapy = as.factor(chemotherapy),
              er_status_measured_by_ihc = as.factor(er_status_measured_by_ihc),
              her2_status_measured_by_snp6 = as.factor(her2_status_measured_by_snp6),
              hormone_therapy = as.factor(hormone_therapy),
              inferred_menopausal_state = as.factor(inferred_menopausal_state),
              integrative_cluster = as.factor(integrative_cluster),
              pam50_claudin_low_subtype = as.factor(pam50_claudin_low_subtype),
              three_gene_classifier_subtype = as.factor(three_gene_classifier_subtype),
              patients_vital_status = as.factor(patients_vital_status),
              primary_tumor_laterality = as.factor(primary_tumor_laterality),
              radio_therapy = as.factor(radio_therapy),
              tumor_other_histologic_subtype = as.factor(tumor_other_histologic_subtype),
              type_of_breast_surgery = as.factor(type_of_breast_surgery),
              # numeric features
              age_at_diagnosis = as.numeric(age_at_diagnosis), 
              lymph_nodes_examined_positive = as.numeric(lymph_nodes_examined_positive),
              nottingham_prognostic_index = as.numeric(nottingham_prognostic_index),
              cohort = as.numeric(cohort),
              overall_survival_months = as.numeric(overall_survival_months), 
              relapse_free_status_months = as.numeric(relapse_free_status_months) 
              ) #%>% 


# loop to check and create masks for all columns in cp with missing data
for (col in colnames(cp)) {
  if (any(is.na(cp[[col]]))) {
    mask_col <- paste0(col, "_mask")
    cp[[mask_col]] <- if_else(!is.na(cp[[col]]), 1, 0)
  }
}

sapply(cp, function(x) sum(is.na(x)) )

# merge clinical sample and patient data
dt <- merge(cs, cp, by = "patient_id", all.x = TRUE)
names(dt)
 
#################################################################################
# mRNA expression data ----
#################################################################################

# mexp: first column is gene symbol, second column is entrez id, rest are sample IDs 
# colnames are [1] "Hugo_Symbol"    "Entrez_Gene_Id" "MB-0362"        "MB-0346" etc
# not z-scores, likely log2 transformed intensities 
genes <- mexp$Hugo_Symbol
exp_mat <- as.matrix(mexp[ , -c(1,2), with = FALSE])

# simply drop genes with > 20% missing values
na_by_gene <- rowMeans(is.na(exp_mat))
exp_mat <- exp_mat[na_by_gene <= 0.2, , drop=FALSE]

# simply impute remaining NAs gene-wise by median
if (anyNA(exp_mat)) {
  meds <- rowMedians(exp_mat, na.rm = TRUE)
  idx <- which(is.na(exp_mat), arr.ind = TRUE)
  exp_mat[idx] <- meds[idx[,1]]
}

rownames(exp_mat) <- genes
colnames(exp_mat) <- colnames(mexp)[-c(1,2)]

# common samples between clinical and expression data and filter to only those samples
common_samples <- intersect(colnames(exp_mat), dt$patient_id) # 1980 samples
exp_mat <- exp_mat[ , common_samples, drop = FALSE]


dt2 <- dt %>%
  dplyr::filter(patient_id %in% common_samples) %>%
  dplyr::arrange(match(patient_id, colnames(exp_mat)))

grp <- factor(dt2$tumor_stage, levels = c("early", "late")) # early = 1338; late = 128; n = 1466 samples
table(grp)

# check whether exp_mat are already normalized log2 (range is 0 to 16)
summary(as.numeric(exp_mat))

# Pre-processing steps for mRNA expression data
## Train, val, test split before preprocessing
set.seed(202501)                               # reproducibility
sample_ids <- colnames(exp_mat)
n_total <- length(sample_ids)
train_frac <- 0.70
val_frac <- 0.15

train_ids <- sample(sample_ids, size = floor(train_frac * n_total))
rem_ids <- setdiff(sample_ids, train_ids)
val_ids <- sample(rem_ids, size = floor(val_frac * n_total))
test_ids <- setdiff(rem_ids, val_ids)

train_mat_raw <- exp_mat[, train_ids, drop = FALSE]  
val_mat_raw <- exp_mat[, val_ids, drop = FALSE]
test_mat_raw <- exp_mat[, test_ids, drop = FALSE]



## Preprocess
qn_reference <- compute_qn_reference(train_mat_raw)
train_qn <- apply_qn(train_mat_raw, qn_reference)
val_qn <- apply_qn(val_mat_raw, qn_reference)
test_qn <- apply_qn(test_mat_raw, qn_reference)

train_log2 <- safe_log2(train_qn)
val_log2 <- safe_log2(val_qn)
test_log2 <- safe_log2(test_qn)

gene_means <- rowMeans(train_log2, na.rm = TRUE)
gene_sds <- apply(train_log2, 1, sd, na.rm = TRUE)
gene_sds[is.na(gene_sds) | gene_sds == 0] <- 1  # guard against zero variance

#expr_train <- zscore_with_stats(train_log2, gene_means, gene_sds)
#expr_val <- zscore_with_stats(val_log2, gene_means, gene_sds)
#expr_test <- zscore_with_stats(test_log2, gene_means, gene_sds)


expr_train <- train_log2
expr_val <- val_log2
expr_test <- test_log2


expr_splits <- list(
  train = expr_train,
  val = expr_val,
  test = expr_test,
  split_ids = list(train = train_ids, val = val_ids, test = test_ids),
  stats = list(qn_reference = qn_reference, means = gene_means, sds = gene_sds)
)

clinical_splits <- lapply(expr_splits$split_ids, function(ids) align_clinical_split(dt2, ids))
stopifnot(
  identical(colnames(expr_splits$train), clinical_splits$train$patient_id),
  identical(colnames(expr_splits$val), clinical_splits$val$patient_id),
  identical(colnames(expr_splits$test), clinical_splits$test$patient_id)
)

## Diagnostic plots on train set
train_plot_df <- data.frame(
  value = c(as.vector(train_mat_raw), as.vector(expr_train)),
  state = rep(c("Raw", "Processed"),
              times = c(length(train_mat_raw), length(expr_train)))
) %>%
  dplyr::filter(!is.na(value))

p_train_raw <- train_plot_df %>%
  dplyr::filter(state == "Raw") %>%
  ggplot2::ggplot(ggplot2::aes(x = value)) +
  ggplot2::geom_density(fill = "#F8766D", alpha = 0.6) +
  ggplot2::labs(title = "Train distribution before preprocessing",
                x = "Expression", y = "Density") +
  ggplot2::theme_minimal()

p_train_processed <- train_plot_df %>%
  dplyr::filter(state != "Raw") %>%
  ggplot2::ggplot(ggplot2::aes(x = value)) +
  ggplot2::geom_density(fill = "#00BFC4", alpha = 0.6) +
  ggplot2::labs(title = "Train distribution after preprocessing",
                x = "Z-score", y = "Density") +
  ggplot2::theme_minimal()

train_preprocessing_check <- cowplot::plot_grid(p_train_raw, p_train_processed, ncol = 2)

#################################################################################
## DEG identification using limma-voom and edgeR on training set only
###################################################################################
## keep genes with at least 20% samples having expression > 1 (z-score)
keep_genes <- apply(expr_train, 1, function(x) sum(x > 1) >= 0.2 * length(x))
dt_train_all <- clinical_splits$train
stage_complete <- !is.na(dt_train_all$tumor_stage)
dt_train <- dt_train_all[stage_complete, , drop = FALSE]

expr_train_stage <- expr_train[, dt_train$patient_id, drop = FALSE]
expr_filt <- expr_train_stage[keep_genes, , drop = FALSE]

grp_train <- factor(dt_train$tumor_stage, levels = c("early", "late"))
stopifnot(length(grp_train) == ncol(expr_filt))

design <- stats::model.matrix(~ grp_train) # keep groups for IDs in train set only
fit <- limma::lmFit(expr_filt, design)
fit <- limma::eBayes(fit)
fit$genes <- data.frame(gene = rownames(expr_filt))

## limma results: late vs early [add p-value adjustment p.adjust, method = "BH"]
res_limma <- topTable(fit, coef = "grp_trainlate", number = Inf) %>%
  as.data.frame()

## select top up and down regulated genes in late stage breast cancer
deg_limma <- res_limma %>%
  dplyr::mutate(direction = ifelse(logFC > 0, "up", "down"))

deg_limma_selected <- deg_limma %>%
  dplyr::filter(adj.P.Val < 0.05, abs(logFC) >= 1) %>%
  dplyr::arrange(adj.P.Val)

if (nrow(deg_limma_selected) == 0) {
  warning("No DEGs passed adj.P.Val < 0.05 and |logFC| >= 1; using top 50 genes by adjusted p-value instead.")
  deg_limma_selected <- deg_limma %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = 50)
}

deg_genes_limma <- deg_limma_selected$gene
deg_genes <- deg_genes_limma  
length(deg_genes)

# create latex table of DEGs from limma results (top 10)
deg_tbl_up <- deg_limma_selected %>%
  dplyr::filter(direction == "up") %>%
  dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
  head(10)

deg_tbl_down <- deg_limma_selected %>%
  dplyr::filter(direction == "down") %>%
  dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
  tail(10)

deg_tbl <- rbind(deg_tbl_up, deg_tbl_down)
knitr::kable(deg_tbl, format = "latex", booktabs = TRUE, digits = 4,
      caption = "Top 10 Upregulated and Downregulated DEGs identified by limma-voom method")

# Volcano plot for limma results [improve this later]
volcano_data <- res_limma %>%
  dplyr::mutate(
    significant = case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE                          ~ "Not significant"
    ),
    neglog10_adjP = -log10(adj.P.Val)
  )

## genes to label for every class
## Up/Down: label top N by abs logFC
## Not significant: label top N by "volcano extremeness" (high -log10(p) AND large abs(logFC))
top_n_sig  <- 10
top_n_nsig <- 10

label_data <- bind_rows(
  volcano_data %>%
    dplyr::filter(significant %in% c("Upregulated", "Downregulated")) %>%
    dplyr::arrange(desc(abs(logFC))) %>%
    dplyr::group_by(significant) %>%
    dplyr::slice_head(n = top_n_sig) %>%
    dplyr::ungroup(),
  
  volcano_data %>%
    dplyr::filter(significant == "Not significant") %>%
    dplyr::mutate(extremeness = neglog10_adjP * abs(logFC)) %>%  # simple, stable ranking
    dplyr::arrange(desc(extremeness)) %>%
    dplyr::slice_head(n = top_n_nsig)
)

volcano_plot <- ggplot2::ggplot(volcano_data, aes(x = logFC, y = neglog10_adjP, color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "#969696"
  ) +
  scale_color_manual(values = c(
    "Upregulated"     = "red",
    "Downregulated"   = "blue",
    "Not significant" = "gray"
  )) +
  labs(
    title = "Volcano plot of DEGs (limma)",
    x = "Log2 Fold Change (Late vs Early)",
    y = "-Log10 Adjusted P-Value",
    color = "DEG class"
  ) +
  theme_minimal()
volcano_plot

# need to do sensitivity analysis with edgeR as well
top_genes <- c(deg_tbl_up$gene, deg_tbl_down$gene)
top_genes <- top_genes[top_genes %in% rownames(expr_splits$train)]

if (!length(top_genes)) {
  stop("No DEG genes overlap with expression matrix after filtering; cannot proceed to downstream analyses.")
}


# Downstream analysis data

## Expression of top DEGs
exp_top <- expr_splits$train[top_genes, dt_train$patient_id, drop = FALSE]
exp_top_df <- as.data.frame(t(exp_top)) %>%
  dplyr::mutate(patient_id = rownames(.)) %>%
  dplyr::left_join(dt_train %>% dplyr::select(patient_id, tumor_stage), by = "patient_id") %>%
  tidyr::pivot_longer(cols = all_of(top_genes), names_to = "gene", values_to = "expression")

## Boxplots of top DEGs expression by tumor stage 
boxplot_top_degs <- ggplot(exp_top_df, aes(x = tumor_stage, y = expression, fill = tumor_stage)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = "free_y") +
  labs(title = "Expression of Top DEGs by Tumor Stage",
       x = "Tumor Stage",
       y = "Expression (Z-score)") +
  theme_minimal() +
  theme(legend.position = "none")
boxplot_top_degs

#################################################################################
# Analysis methods ----
################################################################################

## (1) DR Penalized Logistic Regression ----

## Fit on train set and test on test set computing AUC for tumor stage prediction then estimate ATEs for each DEG adjusting for 
# confounders in the dataset dt

dr_covariates <- c(
  "age_at_diagnosis",
  "er_status",
  "her2_status",
  "neoplasm_histologic_grade",
  "lymph_nodes_examined_positive",
  "nottingham_prognostic_index"
)

confounder_df <- dt_train %>% dplyr::select(all_of(dr_covariates))
complete_mask <- complete.cases(confounder_df) & !is.na(dt_train$tumor_stage)

if (!any(complete_mask)) {
  stop("No complete cases available for DR estimation after filtering confounders/outcome.")
}

# Restrict to subjects with fully observed confounders so glmnet receives a clean design matrix.
dt_dr <- dt_train[complete_mask, , drop = FALSE]
expr_train_dr <- expr_splits$train[, dt_dr$patient_id, drop = FALSE]
confounder_mm <- model.matrix(~ . - 1, data = confounder_df[complete_mask, , drop = FALSE])
X_dr <- as.data.frame(confounder_mm)
Y_dr <- ifelse(dt_dr$tumor_stage == "late", 1, 0)

## Helper to fit glmnet and predict [move to helpers]
fit_glmnet_cv <- function(X, y, family = "gaussian", alpha = 1) {
  X_mat <- as.matrix(X)
  if (anyNA(X_mat) || anyNA(y)) {
    stop("Input data for glmnet contains NA values after filtering.")
  }
  cvfit <- cv.glmnet(X_mat, y, family = family, alpha = alpha)
  list(
    model = cvfit,
    predict = function(newX) {
      predict(cvfit, newx = as.matrix(newX), s = "lambda.min")[, 1]
    }
  )
}

## DoubleML-style PLR for binary outcome (linear probability model) [move to helpers later]
## Y: outcome vector
## D: treatment vector (gene expression)
## Idea: Fit outcome model Y = theta*D + g(X) + e and treatment model D = m(X) + v using K-fold cross-fitting then regress residuals \tilde{Y} on \tilde{D} to get theta (Chernozhukov et al. (2016)). 
## Here we fit both g(X) and m(X) using Lasso regression via glmnet with cross-validation to select lambda.
dml_plr <- function(Y, D, X, K = 5) {
  n <- length(Y)
  folds <- sample(rep(1:K, length.out = n))
  
  y_tilde <- rep(NA_real_, n)
  d_tilde <- rep(NA_real_, n)
  
  for (k in 1:K) {
    idx_train <- which(folds != k)
    idx_test  <- which(folds == k)
    
    X_train <- X[idx_train, , drop = FALSE]
    X_test  <- X[idx_test,  , drop = FALSE]
    
    Y_train <- Y[idx_train]
    D_train <- D[idx_train]
    
    fit_Y <- fit_glmnet_cv(X_train, Y_train, family = "gaussian")
    g_hat <- fit_Y$predict(X_test)
    
    fit_D <- fit_glmnet_cv(X_train, D_train, family = "gaussian")
    m_hat <- fit_D$predict(X_test)
    
    y_tilde[idx_test] <- Y[idx_test] - g_hat
    d_tilde[idx_test] <- D[idx_test] - m_hat
  }
  
  fit_final <- lm(y_tilde ~ d_tilde)
  theta_hat <- coef(fit_final)["d_tilde"]
  se_theta  <- summary(fit_final)$coefficients["d_tilde", "Std. Error"]
  
  list(
    theta = theta_hat,
    se    = se_theta,
    ci    = c(theta_hat - 1.96 * se_theta,
              theta_hat + 1.96 * se_theta)
  )
}

#  pick top 5 DEGs and run dml_plr on each adjusting for confounders and return gene, theta, se, ci_lower, ci_upper
deg_sub <- top_genes[seq_len(min(5, length(top_genes)))]
if (!length(deg_sub)) {
  stop("No DEGs available for DR estimation step.")
}
results_plr <- list()
for (gene in deg_sub) {
  if (!gene %in% rownames(expr_train_dr)) {
    warning(sprintf("Skipping %s because expression values are unavailable after filtering.", gene))
    next
  }
  gene_expr <- expr_train_dr[gene, , drop = TRUE]
  D <- as.numeric(gene_expr)
  
  dml_result <- dml_plr(Y = Y_dr, D = D, X = X_dr, K = 5)
  
  results_plr[[gene]] <- data.frame(
    gene = gene,
    theta = dml_result$theta,
    se = dml_result$se,
    ci_lower = dml_result$ci[1],
    ci_upper = dml_result$ci[2]
  )
}

results_plr_df <- do.call(rbind, results_plr)
results_plr_df



