library(DESeq2)
cds = DESeqDataSetFromMatrix(countData = count_matrix,
                             colData = expgroup,
                             design = ~cancer+oxygen)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds)
cds = DESeq(cds)
res = results(cds)

diffexpgenes = rownames(subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)))
normvalues = counts(cds, normalized = T)
diffexpvalues = data.frame(normvalues[diffexpgenes,])
clusters = hclust(dist(diffexpvalues))
groups = cutree(clusters, k = 3)
print(table(groups))

library(pheatmap)
pheatmap(diffexpvalues, 
         cluster_rows = TRUE, 
         scale = "row", annotation_col = as.data.frame(expgroup),
         annotation_row = as.data.frame(groups)["groups"])
