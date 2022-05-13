Triple-Negative Breast Cancer (TNBC) is a highly aggressive tumor subtype with a high risk of recurrence, metastasis, chemotherapy resistance, and acquired capacity to survive and grow under nutrient-deprived and hypoxic conditions. Under hypoxic conditions, TNBC cells can grow, survive, induce metabolic reprogramming and apoptosis and alter cell adhesion and motility to facilitate metastasis and resistance to chemotherapy. The project aims to identify what cellular organization genes are upregulated by tumor cells in hypoxic and/or mTOR inhibited environments and if they are a part of the same pathways.

The dataset includes expression profiling by high throughput sequencing. The study consists of only the total RNA and excludes Ribosomal RNA. Thus, the data consists of sixteen samples with two replicates per sample. The samples belong to either the Benign or Malignant (TNBC) group. Each sample was exposed to Normoxic and Hypoxic environments and then further treated with an mTOR pathway inhibitor, PP242. This produced eight different sample groups leading to a multifactorial analysis.
The SRA accession list consisting of 16 paired-end samples was converted to their corresponding 32 FASTQ files. FASTQC analysis was performed to check the quality of the FASTQ files. No trimming was performed since the alignment was 95% accurate with no adapter sequences.

These FASTQ files were aligned using HISAT and converted to corresponding 16 SAM files. SAMTOOLS converted the SAM files to BAM format for efficient storage and easy manipulation. The BAM files were sorted and indexed to aid in the quick extraction of alignments overlapping particular genomic regions. The counts from each BAM file were generated using HTSeq and merged to receive the entire count matrix with 16 samples and 64,256 genes. 

The differential gene expression analysis was performed using DESeq to understand the biological differences between healthy and diseased states under varying conditions. To test the hypothesis, three experimental design groups were created:

1.	Cancer + Oxygen
2.	Cancer + Treatment
3.	Cancer + Oxygen + Treatment

The dispersion estimates were plotted to understand the dispersion with respect to the mean of the normalized count for each experimental design. The differentially expressed genes were identified for each experimental design, and the normvalues and diffexpvalues were calculated. 

Hierarchical clustering was performed to analyze the tree structure from the data similarity, and the relationship between different sub-clusters was observed. Heatmaps were plotted to analyze and visualize this multi-dimensional data to gain a profound understanding of the high-throughput gene expression data.  

The GO Term Enrichment Analysis is performed to classify gene sets and receive a holistic understanding through the corresponding annotations. The ENSEMBL gene ids were converted to ENTREZ ids using biomaRt. The GO term enrichment analysis was performed using the hyperGTest function.
