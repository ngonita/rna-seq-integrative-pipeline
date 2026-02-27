# ==============================================================================
# RNA-SEQ + METHYLATION INTEGRATIVE ANALYSIS
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(123)

library(edgeR)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)

dir.create("results", showWarnings = FALSE)

# -----------------------------
# INPUT VALIDATION
# -----------------------------
if(!file.exists("data/counts_matrix.txt"))
  stop("Counts matrix missing")

if(!file.exists("data/metadata.txt"))
  stop("Metadata missing")

# -----------------------------
# LOAD DATA
# -----------------------------
counts_matrix <- read.table(
  "data/counts_matrix.txt",
  header=TRUE,
  row.names=1
)

metadata <- read.table(
  "data/metadata.txt",
  header=TRUE,
  row.names=1
)

if(!all(colnames(counts_matrix)==rownames(metadata))){
  stop("Mismatch between metadata and counts")
}

# -----------------------------
# RNA-SEQ ANALYSIS
# -----------------------------
dge <- DGEList(counts=counts_matrix,
               group=metadata$Condition)

keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~Condition, data=metadata)

v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

results_RNA <- topTable(
  fit,
  coef=2,
  number=Inf,
  adjust.method="BH"
)

write.csv(results_RNA,
          "results/RNAseq_DE_results.csv")

sig_genes_RNA <- results_RNA %>%
  filter(adj.P.Val < 0.05 &
         abs(logFC) > 1)

# -----------------------------
# METHYLATION DATA
# -----------------------------
results_Methy <- read.table(
  "data/res_methylation.txt",
  header=TRUE
)

# -----------------------------
# INTEGRATION
# -----------------------------
merged_data <- inner_join(
  results_RNA,
  results_Methy,
  by="GeneName",
  suffix=c("_RNA","_Methy")
)

final_signature <- merged_data %>%
  filter(adj.P.Val_RNA < 0.05 &
         adj.P.Val_Methy < 0.05) %>%
  filter(
    (logFC_RNA>0 & logFC_Methy<0) |
    (logFC_RNA<0 & logFC_Methy>0)
  )

write.csv(final_signature,
          "results/final_signature.csv")

# -----------------------------
# VISUALIZATION
# -----------------------------
p <- ggplot(merged_data,
       aes(logFC_Methy,logFC_RNA))+
  geom_point(alpha=0.4,color="grey")+
  geom_point(data=final_signature,
             color="red",size=2)+
  geom_hline(yintercept=0,
             linetype="dashed")+
  geom_vline(xintercept=0,
             linetype="dashed")+
  theme_minimal()

ggsave("results/starburst_plot.png",
       p,width=6,height=5)

# -----------------------------
# SESSION TRACEABILITY
# -----------------------------
writeLines(
  capture.output(sessionInfo()),
  "results/sessionInfo.txt"
)
