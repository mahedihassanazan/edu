library(ggplot2)
library(ggrepel)

cat("Reading data...\n")
res <- read.csv("Retina/Retina_DEGs_Result.csv", row.names=1)

# Handle NA values in padj (DESeq2 sets some to NA if baseMean is too low)
res$padj[is.na(res$padj)] <- 1

# Create a new column for grouping colors
res$Expression <- "Not Significant"
res$Expression[res$padj < 0.05 & res$log2FoldChange > 1] <- "Upregulated"
res$Expression[res$padj < 0.05 & res$log2FoldChange < -1] <- "Downregulated"

# Extract just the gene symbol from rownames (e.g. "ENSMUSG00000025902|Sox17" -> "Sox17")
res$Symbol <- sapply(strsplit(rownames(res), "\\|"), function(x) x[2])

# Find the top 10 Upregulated and top 10 Downregulated genes by p-value to label them
top_up <- head(res[res$Expression == "Upregulated", ][order(res[res$Expression == "Upregulated", ]$padj), ], 10)
top_down <- head(res[res$Expression == "Downregulated", ][order(res[res$Expression == "Downregulated", ]$padj), ], 10)
top_genes <- rbind(top_up, top_down)

cat("Generating Volcano Plot...\n")
p <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=Expression)) +
    geom_point(alpha=0.7, size=1.5) +
    scale_color_manual(values=c("Downregulated"="#2b8cbe", "Not Significant"="#bdbdbd", "Upregulated"="#e34a33")) +
    geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black", alpha=0.5) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", alpha=0.5) +
    geom_text_repel(data=top_genes, aes(label=Symbol), size=4, color="black", max.overlaps=Inf) +
    theme_minimal() +
    coord_cartesian(xlim = c(-12, 12)) +
    labs(title="Volcano Plot of DEGs (Glaucoma vs Normal)",
         x="log2 Fold Change",
         y="-log10(Adjusted P-value)") +
    theme(legend.position="right",
          plot.title = element_text(hjust = 0.5, size=16, face="bold"),
          axis.title = element_text(size=14),
          legend.title = element_blank(),
          legend.text = element_text(size=12))

png("Volcano_Plot.png", width=1200, height=1000, res=150)
print(p)
dev.off()

cat("Done! Plot saved as Volcano_Plot.png\n")
