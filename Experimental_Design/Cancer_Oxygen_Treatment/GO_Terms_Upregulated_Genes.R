library(DESeq2)
cds = DESeqDataSetFromMatrix(countData = count_matrix,
                             colData = expgroup,
                             design = ~cancer+oxygen+treatment)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds)
cds = DESeq(cds)
res = results(cds)
diffexpgenes = rownames(subset(res, padj < 0.05 & (log2FoldChange > 1)))
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

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
diffexp_genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=diffexpgenes,
  mart=mart)

params = new("GOHyperGParams", 
             geneIds = diffexp_genes$entrezgene_id, 
             universeGeneIds = genes$entrezgene_id,
             annotation = "org.Hs.eg.db",
             ontology = "BP",
             pvalueCutoff = 0.001,
             testDirection = "over")

test = hyperGTest(params)
print(summary(test))

test_df = data.frame(summary(test))
save_df = data.frame(test_df$Term, test_df$Pvalue, test_df$Count, test_df$Size)
write.table(save_df, "COT_U_GOTerms.csv", sep=",")
