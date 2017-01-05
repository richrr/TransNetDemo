#library(TransNetDemo)
library(stringr)
library(ProNet)
library(igraph)


individualPvalueCutoff = 0.3
combinedPvalueCutoff = 0.05
combinedFDRCutoff = 0.2

# groups to compare
groupA = "HFHS"
groupB = "NCD"
sampleIdColName = "SampleID"
factorColName = "Factor"
geneSymbolColName ="IdSymbol"
# fold change column based on mean or median
foldchVar = "FoldChange_HFHS_NCD"   # or "FoldChange_Median_HFHS_NCD"

# map files
mapf1 = system.file("extdata", "mapping_file.rand.1.tsv", package = "TransNetDemo")
mapf2 = system.file("extdata", "mapping_file.rand.2.tsv", package = "TransNetDemo")

# gene files
genef1 = system.file("extdata", "gene_file_1.tsv", package = "TransNetDemo")
genef2 = system.file("extdata", "gene_file_2.tsv", package = "TransNetDemo")

# genes
Comp_genes = Compare_groups(mapf1, genef1, sampleIdColName, factorColName , groupA, groupB, geneSymbolColName)

# Select differentially abundant elements:
#########  (i) using a significance threshold #########
Sign_genes = Comp_genes[which(Comp_genes$FDR < 0.05),]
# this data has been saved as "Sign_genes_precomputed_1"
######### ######### ######### #########
# OR
######### (ii) in case you have multiple datasets, we recommend meta-analysis ###########
## we have already run these steps and stored the data under the variable "Sign_genes_metaanalysis_precomputed" ##
Comp_genes1 = Comp_genes
Comp_genes2 = Compare_groups(mapf2, genef2, sampleIdColName, factorColName , groupA, groupB, geneSymbolColName)
numbDatasets=2  # number of datasets

s_df = merge(Comp_genes1, Comp_genes2, by="row.names")
rownames(s_df) = rownames(Comp_genes1)
s_df = Check_consistency(s_df, foldchVar, 1, numbDatasets)
comb_in_df = Calc_combined(s_df)
Sign_genes = Apply_sign_cutoffs(comb_in_df, individualPvalueCutoff, combinedPvalueCutoff, combinedFDRCutoff )

# if this returns TRUE, you did everything correct!
identical(rownames(Sign_genes), rownames(Sign_genes_metaanalysis_precomputed))
#write.csv(Sign_genes,"Sign_genes_File.csv", quote=FALSE)
######### ######### ######### #########


# gene pairs
Corr_pairs = Correlation_in_group(mapf1, genef1, sampleIdColName, factorColName, groupA,  geneSymbolColName, rownames(Sign_genes))

# Select significant correlations:
#########  (i) using a significance threshold #########
Sign_pairs = Corr_pairs[which(Corr_pairs$pvalue < 0.05 & Corr_pairs$FDR < 0.1),]
######### ######### ######### #########
# OR
######### (ii) in case you have multiple datasets, we recommend meta-analysis ###########
## we have already run these steps and stored the data under the variable "Sign_genepairs_metaanalysis_precomputed" ##
Corr_pairs1 = Corr_pairs
Corr_pairs2 = Correlation_in_group(mapf2, genef2, sampleIdColName, factorColName, groupA,  geneSymbolColName, rownames(Sign_genes))

s_df = merge(Corr_pairs1, Corr_pairs2, by="row.names")
rownames(s_df) = rownames(Corr_pairs1)
s_df = Check_consistency(s_df, "Coefficient", 0, numbDatasets)
comb_in_df = Calc_combined(s_df)
Sign_pairs = Apply_sign_cutoffs(comb_in_df,individualPvalueCutoff, combinedPvalueCutoff, combinedFDRCutoff )

# if this returns TRUE, you did everything correct!
identical(rownames(Sign_pairs), rownames(Sign_genepairs_metaanalysis_precomputed))
#write.csv (Sign_pairs,"Sign_genes_pairs_File.csv", quote=FALSE)
######### ######### ######### #########


genes_df = Calc_median_val(Sign_genes, foldchVar)
identical(genes_df, Gene_df_precomputed)
pairs_df = Calc_median_val(Sign_pairs, "Coefficient")
outNetwork = Puc_compatible_network(pairs_df, genes_df)

# if this returns TRUE, you did everything correct!
identical(outNetwork[,c("partner1", "partner2")], Gene_network_precomputed[,c("partner1", "partner2")])
#write.csv (outNetwork,"gene-networkFile.csv", quote=FALSE)

cluster1 = Identify_subnetworks(outNetwork)
summary(cluster1)

# if this returns TRUE, you did everything correct!
identical(get.edgelist(cluster1), get.edgelist(genes_mcode_cluster1_precomputed)) # OR #sum(get.adjacency(cluster1) != get.adjacency(genes_mcode_cluster1_precomputed)) # should return 0
#write_graph(cluster1, "genes_mcode_cluster1_edges.txt", "ncol")

plot(cluster1, vertex.label=NA, layout= layout_in_circle, vertex.size=3)
