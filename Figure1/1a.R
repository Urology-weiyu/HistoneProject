# ============================================================
# Title:    Pan-Cancer Correlation Analysis of LDHA Expression 
#           and Cuproptosis-Inhibitory Signature
# Description: 
# Perform a pancancer-wide Pearson correlation analysis between LDHA expression 
# and the ssGSEA-derived "Cuproptosis inhibitory" gene set activity across TCGA datasets. 
# The script integrates expression, clinical, and gene set enrichment data, and 
# visualizes correlation coefficients across tumor types as a heatmap.
# ============================================================

# ---------- [1. Environment Setup] ----------
setwd('path/to/workdir')
library(dplyr)
library(purrr)
library(GSEABase)
library(clusterProfiler)
library(GSVA)
library(ggplot2)

# ---------- [2. Load pancancer data] ----------
load('path/to/TCGA_pancancer_clin.rdata')
load('path/to/mRNAexp_TPM_pancancer_33.rdata')

# ---------- [3. Extract expression for target gene] ----------
inputgene <- c('LDHA')
exp <- expall[inputgene, ] %>%
  t() %>%
  as.data.frame()
exp$sample_id <- substr(rownames(exp), 1, 15)

# Merge with TCGA clinical metadata
tmp <- merge(exp, tcga_clin, by = 'sample_id')
tmp$sample_type <- substr(tmp$sample_id, 14, 15)
tmp <- subset(tmp, sample_type == '01')  # tumor only
tmp <- subset(tmp, select = c(1:(length(inputgene) + 3)))

# ---------- [4. Compute ssGSEA score for Cuproptosis-Inhibitory signature] ----------
geneset <- read.gmt('path/to/Cuinhib.gmt')
es.max <- gsva(as.matrix(expall),
               geneset,
               method = 'ssgsea',
               mx.diff = FALSE,
               verbose = FALSE,
               parallel.sz = 1)

group <- es.max %>%
  t() %>%
  as.data.frame()
group$id <- rownames(group)
colnames(group) <- c('Cuproptosis_inhibitory', 'sample_id')
group$id <- substr(group$sample_id, 1, 15)

# ---------- [5. Merge expression and ssGSEA results] ----------
df <- merge(tmp, group, by = 'sample_id')

# ---------- [6. Compute Pearson correlation per cancer type] ----------
cor1 <- df %>%
  group_by(project) %>%
  group_split() %>%
  map_df(~{
    test <- cor.test(.x$Cuproptosis_inhibitory, log10(.x$LDHA), method = "pearson")
    tibble(
      project = unique(.x$project),
      correlation = test$estimate,
      p_value = test$p.value
    )
  })

cor1$project <- factor(cor1$project,
                       levels = rev(cor1$project[order(cor1$correlation)]))

# ---------- [7. Plot correlation heatmap] ----------
p <- ggplot(cor1, aes(
  y = "LDHA-\nCuproptosis\nInhibitory",
  x = project,
  fill = correlation
)) +
  geom_tile(color = "black", width = 0.9, height = 0.9, size = 0.5) +
  scale_fill_gradient2(
    low = "#1a237e",
    mid = "white",
    high = "#b71c1c",
    midpoint = 0
  ) +
  labs(x = NULL, y = '', title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 13, color = 'black', hjust = .5),
    axis.text.x = element_text(size = 13, color = 'black'),
    axis.title.y = element_text(size = 14),
    axis.ticks.length = unit(2, "mm"),
    legend.position = 'none'
  ) +
  coord_fixed()

ggsave('Fig1a.pdf',
       plot = p,
       width = 8,
       height = 3.5)
