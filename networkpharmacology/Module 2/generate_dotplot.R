options(repos = c(CRAN = "https://cloud.r-project.org"))
options(scipen=999) # This forces standard decimal format everywhere!

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

cat("Reading data...\n")
res <- read.csv("Retina/Retina_DEGs_Result.csv", row.names=1)

# ----- UPREGULATED -----
cat("Filtering Upregulated DEGs...\n")
res_up <- res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, ]
symbols_up <- sapply(strsplit(rownames(res_up), "\\|"), function(x) x[2])
go_up <- enrichGO(gene = symbols_up, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

cat("Generating Upregulated Plot...\n")
png("GO_DotPlot_Upregulated.png", width=1200, height=1200, res=150)
print(dotplot(go_up, showCategory=10, label_format=40) + 
      ggtitle("Top 10 Biological Processes (Upregulated)"))
dev.off()

# ----- DOWNREGULATED -----
cat("Filtering Downregulated DEGs...\n")
res_down <- res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, ]
symbols_down <- sapply(strsplit(rownames(res_down), "\\|"), function(x) x[2])
go_down <- enrichGO(gene = symbols_down, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

cat("Generating Downregulated Plot...\n")
png("GO_DotPlot_Downregulated.png", width=1200, height=1200, res=150)
print(dotplot(go_down, showCategory=10, label_format=40) + 
      ggtitle("Top 10 Biological Processes (Downregulated)"))
dev.off()

cat("All done!\n")
