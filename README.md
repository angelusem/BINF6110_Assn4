# BINF6110 Assignment 4: scRNA-seq analysis of influenza-restricted nasal mucosa

## Introduction

In this analysis, a Seurat object derived from a mouse influenza infection study of the nasal mucosa was analyzed. The dataset includes multiple nasal compartments, multiple timepoints after infection, and mouse-level metadata, making it suitable for clustering, annotation, and cluster-specific differential expression analysis.

The main goals of this analysis were to identify major cell populations in the dataset, annotate clusters using both automated and marker-based approaches, and evaluate transcriptional changes within at least one cluster across experimental groups. The primary downstream comparison focused on a macrophage-supported cluster in respiratory mucosa (RM), comparing early infection (D02) with baseline (Naive). A second, more cleanly annotated endothelial cluster was also analyzed as a supporting comparison.

## Methods

### Data input and project setup

The analysis was performed in RStudio within an `renv`-managed project environment.

The main processed object generated during analysis was:

-   `data/seurat_ass4_lognorm_clustered.rds`

### Quality control and preprocessing

The initial object was inspected for metadata structure, tissue distribution, and timepoint distribution. Mitochondrial percentage was calculated using the mouse-appropriate `^mt-` prefix. Quality control was assessed using violin plots and scatter plots of `nFeature_RNA`, `nCount_RNA`, and `percent.mt`. Cells with missing `mouse_id`, fewer than 200 detected features, or mitochondrial fraction greater than 20% were removed.

A log-normalization workflow was employed, consisting of `NormalizeData()`, `FindVariableFeatures()`, `ScaleData()` with `percent.mt` regressed out, `RunPCA()`, `FindNeighbors()`, `FindClusters()`, and `RunUMAP()`. Clustering yielded 36 clusters. Harmony or another integration method was considered but not applied, because the UMAP showed biologically structured clustering rather than strong sample-driven separation.

### Cluster annotation

Cluster annotation was performed using two complementary approaches. First, Seurat marker detection was used to identify cluster-enriched genes. Second, SingleR was run at the cluster level using `MouseRNAseqData()` as the reference. Marker-based feature plots were then used to validate broad biological compartments, including neuronal, endothelial, myeloid/macrophage, stromal/fibroblast, and epithelial populations. The final labels were kept broad, since this level of annotation was best supported by both the marker evidence and the reference-based labels.

### Differential expression design

For downstream analysis, **cluster 1** was selected as the primary case study because both marker genes and SingleR supported a macrophage identity, and the RM Naive versus D02 comparison provided balanced mouse-level replication with adequate cell numbers per mouse. As a supporting comparison, **cluster 4** was analyzed because it showed a cleaner endothelial identity supported by both marker genes and SingleR.

Cluster-specific differential expression was performed using a **pseudobulk** workflow. Cells were subset by cluster, tissue, and timepoint, then aggregated by `mouse_id` using `AggregateExpression()`. Mouse identity was used as the biological replicate unit. Differential expression was then performed with DESeq2 using the design `~ time`, with the main contrast defined as **D02 vs Naive**, so positive log2 fold change indicates higher expression in D02 and negative log2 fold change indicates higher expression in Naive. Log2 fold changes were then shrunk using `apeglm`, and MA and volcano plots were generated. Because many top-ranked genes in the full DE tables were ribosomal or mitochondrial, filtered interpretation tables excluding genes beginning with `Rpl`, `Rps`, and `mt-` were generated in addition to the full official DE tables.

### Functional analysis

Functional interpretation was performed using ORA and GSEA with `clusterProfiler` and `org.Mm.eg.db`. ORA was run on filtered sets of upregulated genes, while GSEA was run on the full ranked shrunk log2 fold change list after gene symbol to Entrez ID mapping. This allowed both threshold-based and ranked-list-based views of cluster-specific transcriptional change.

## Results

### Overall clustering and annotation

After filtering and log-normalization-based clustering, the dataset separated into 36 clusters with clear broad biological structure. Marker-based visualization and SingleR supported major compartments including neurons, macrophages, endothelial cells, fibroblasts, epithelial cells, B cells, T cells, NK cells, monocytes, and granulocytes. Because the dataset and reference-based annotations most strongly supported broad identities rather than fine subtype mapping, the annotation strategy remained conservative throughout the analysis.

-   **Figure 1.** UMAP colored by cluster or broad SingleR label.\
-   **Figure 2.** Marker validation plot (dot plot or feature plot panel) showing representative neuronal, endothelial, macrophage/myeloid, stromal, and epithelial markers.

### Main downstream comparison: cluster 1 macrophages in RM, D02 vs Naive

Cluster 1 was selected as the main downstream comparison because marker genes and SingleR supported a macrophage identity, and the RM Naive versus D02 contrast provided balanced mouse-level replication. Within RM, the cluster 1 subset contained three biological replicates per group, with D02 counts of 449, 193, and 208 cells and Naive counts of 481, 321, and 264 cells.

Pseudobulk DESeq2 analysis identified a clear transcriptional shift between Naive and D02. In this comparison:

-   **932 genes** had `padj < 0.05`
-   **131 genes** had `padj < 0.05` and `|log2FC| > 1`
-   **104 genes** were up in D02
-   **27 genes** were up in Naive

Many of the most statistically extreme genes in the full table were ribosomal or mitochondrial, so interpretation focused on a filtered set of more biologically interpretable genes. Among the D02-up genes, notable examples included **Plek**, **Ap2a2**, and **Ap2b1**.

ORA on the filtered D02-up gene set returned little support for broad GO Biological Process or Molecular Function themes, but GO Cellular Component identified a narrow repeated signal centered on **clathrin adaptor complex** and **clathrin vesicle coat**, driven mainly by **Ap2a2**, **Ap2b1**, and **Ap1s2**. GSEA on the ranked shrunk log2 fold change list gave a broader signal. On the Naive-up side, enriched terms included **translation-associated processes**, **oxidative phosphorylation**, **respiratory electron transport chain**, and **antigen processing and presentation**. On the D02-up side, enriched terms included **ion transport**, **carbohydrate metabolic process**, **extracellular matrix assembly**, and **export across plasma membrane**.

-   **Figure 3.** MA plot for cluster 1 RM macrophages, D02 vs Naive.\
-   **Figure 4.** Volcano plot for cluster 1 RM macrophages, D02 vs Naive.\
-   **Figure 5.** ORA or GSEA plot for cluster 1 macrophages.\
-   **Figure 6.** Feature plots for selected genes of interest

### Supporting comparison: cluster 4 endothelial cells in RM, D02 vs Naive

Cluster 4 was analyzed as a supporting comparison because it showed a cleaner endothelial identity supported by both SingleR and endothelial marker genes. The RM Naive versus D02 pseudobulk comparison also produced a valid differential expression result:

-   **637 genes** had `padj < 0.05`
-   **111 genes** had `padj < 0.05` and `|log2FC| > 1`
-   **107 interpretable strong-effect genes** remained after removal of ribosomal and mitochondrial genes
-   **82 strong-effect genes** were up in D02
-   **29 strong-effect genes** were up in Naive

Interpretable D02-up genes included **Hes1**, **Cd200**, **Nmt1**. However, the enrichment narrative seemed narrower than in cluster 1. GO Biological Process and Molecular Function ORA returned no significant terms, while GO Cellular Component identified only a small AP-2/clathrin-associated vesicle signal once more driven by **Ap2a2** and **Ap2b1**. GSEA also returned a more limited result, with **mitochondrial translation** enriched on the Naive-up side.

These results made cluster 4 useful as a cleaner comparison cluster, but biologically less rich than cluster 1.

-   **Supplementary Figure 1.** MA or volcano plot for cluster 4 RM endothelial cells, D02 vs Naive.\
-   **Supplementary Figure 2.** Feature plots showing endothelial markers and selected cluster 4 DE genes.

### Overall analytical outcome

Taken together, the analysis recovered broad biologically coherent cell populations and demonstrated that early influenza-associated transcriptional changes could be detected in a cluster-specific manner using mouse-level pseudobulk differential expression. Cluster 1 macrophages provided the richer biological case study, while cluster 4 endothelial cells served as a cleaner but narrower supporting comparison.

## Discussion

The primary conclusion of this analysis is that cluster-specific pseudobulk differential expression recovered meaningful transcriptional changes associated with early infection while preserving biological replication at the mouse level. Cluster 1 macrophages in RM showed the strongest overall case study: the differential expression signal was substantial, the filtered gene list contained interpretable genes, and the enrichment analyses supported changes related to membrane trafficking, transport, extracellular-matrix-associated processes, and mitochondrial or translation-linked programs. Although many top-ranked genes in the full result table were ribosomal or mitochondrial, the combined DE, ORA, and GSEA results support the interpretation that this cluster underwent a real transcriptional state shift between Naive and D02

Cluster 4 endothelial cells helped contextualize the macrophage result. This cluster was easier to annotate, but its enrichment output was narrower and less biologically rich. This contrast was informative because it showed that a cleaner cluster identity did not necessarily produce a stronger pathway-level signal. For this reason, cluster 1 remained the more useful main case study for the assignment, while cluster 4 strengthens the analysis.

Several limitations should also be noted. First, the final cluster annotations were intentionally broad because the dataset and reference-based methods did not justify fine subtype claims with high confidence. Second, local memory limitations prevented SCTransform from being used in the final pipeline, so a log-normalization workflow was adopted instead. Third, many highly ranked DE genes were not ideal narrative anchors, making filtered interpretation tables and pathway-level summaries more informative than raw ranking alone. Finally, enrichment results were uneven across clusters, so the biological interpretation relies on the overall pattern across DE, ORA, and GSEA rather than any single result table.

## References
