library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)

# Load DEGs
df <- read.csv("Simplified_Retina_DEGs_No_NA.csv")

# CRITICAL FIX: Filter for SIGNIFICANT DEGs ONLY (P < 0.05 and |LogFC| > 1)
sig_df <- df %>% filter(P_Value < 0.05 & abs(Log_Fold_Change) > 1)
genes <- as.character(sig_df$Gene_Name)

# Convert Gene Symbol to Entrez ID
gene_ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

# Run KEGG Enrichment
print(paste("Running KEGG Enrichment Analysis on", length(genes), "significant DEGs..."))
kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(kegg_result), "KEGG_Enrichment_Results.csv", row.names=FALSE)

# Generate Barplot
p1 <- barplot(kegg_result, showCategory=15, title="KEGG Pathway Enrichment (Bar Plot)")
ggsave("KEGG_Enrichment_Barplot.png", plot=p1, width=8, height=6, dpi=300)

# Generate Bubble Plot
p2 <- dotplot(kegg_result, showCategory=15, title="KEGG Pathway Enrichment (Bubble Plot)")
ggsave("KEGG_Enrichment_Bubbleplot.png", plot=p2, width=8, height=6, dpi=300)

print("Graphs generated successfully!")
