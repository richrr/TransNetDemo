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
microbeSymbolColName ="IdSymbol"
# fold change column based on mean or median
foldchVar = "FoldChange_HFHS_NCD"   # or "FoldChange_Median_HFHS_NCD"

# map files
mapf1 = system.file("extdata", "mapping_file.rand.1.tsv", package = "TransNetDemo")
mapf2 = system.file("extdata", "mapping_file.rand.2.tsv", package = "TransNetDemo")

# microbe files
microbef1 = system.file("extdata", "microbe_file_1.tsv", package = "TransNetDemo")
microbef2 = system.file("extdata", "microbe_file_2.tsv", package = "TransNetDemo")

# microbes
Comp_microbes = Compare_groups(mapf1, microbef1, sampleIdColName, factorColName , groupA, groupB, microbeSymbolColName)

# Select differentially abundant elements:
#########  (i) using a significance threshold #########
Sign_microbes = Comp_microbes[which(Comp_microbes$FDR < 0.05),]
# this data has been saved as "Sign_microbes_precomputed_1"
######### ######### ######### #########
# OR
######### (ii) in case you have multiple datasets, we recommend meta-analysis ###########
## we have already run these steps and stored the data under the variable "Sign_microbes_metaanalysis_precomputed" ##
Comp_microbes1 = Comp_microbes
Comp_microbes2 = Compare_groups(mapf2, microbef2, sampleIdColName, factorColName , groupA, groupB, microbeSymbolColName)
numbDatasets=2  # number of datasets

s_df = merge(Comp_microbes1, Comp_microbes2, by="row.names")
rownames(s_df) = rownames(Comp_microbes1)
s_df = Check_consistency(s_df, foldchVar, 1, numbDatasets)
comb_in_df = Calc_combined(s_df)
Sign_microbes = Apply_sign_cutoffs(comb_in_df, individualPvalueCutoff, combinedPvalueCutoff, combinedFDRCutoff )

# if this returns TRUE, you did everything correct!
identical(rownames(Sign_microbes), rownames(Sign_microbes_metaanalysis_precomputed))
#write.csv(Sign_microbes,"Sign_microbes_File.csv", quote=FALSE)
######### ######### ######### #########


# microbe pairs
Corr_pairs = Correlation_in_group(mapf1, microbef1, sampleIdColName, factorColName, groupA,  microbeSymbolColName, rownames(Sign_microbes))

# Select significant correlations:
#########  (i) using a significance threshold #########
Sign_pairs = Corr_pairs[which(Corr_pairs$pvalue < 0.05 & Corr_pairs$FDR < 0.1),]
######### ######### ######### #########
# OR
######### (ii) in case you have multiple datasets, we recommend meta-analysis ###########
## we have already run these steps and stored the data under the variable "Sign_microbepairs_metaanalysis_precomputed" ##
Corr_pairs1 = Corr_pairs
Corr_pairs2 = Correlation_in_group(mapf2, microbef2, sampleIdColName, factorColName, groupA,  microbeSymbolColName, rownames(Sign_microbes))

s_df = merge(Corr_pairs1, Corr_pairs2, by="row.names")
rownames(s_df) = rownames(Corr_pairs1)
s_df = Check_consistency(s_df, "Coefficient", 0, numbDatasets)
comb_in_df = Calc_combined(s_df)
Sign_pairs = Apply_sign_cutoffs(comb_in_df,individualPvalueCutoff, combinedPvalueCutoff, combinedFDRCutoff )

# if this returns TRUE, you did everything correct!
identical(rownames(Sign_pairs), rownames(Sign_microbepairs_metaanalysis_precomputed))
#write.csv (Sign_pairs,"Sign_microbes_pairs_File.csv", quote=FALSE)
######### ######### ######### #########


microbes_df = Calc_median_val(Sign_microbes, foldchVar)
identical(microbes_df, Microbe_df_precomputed)
pairs_df = Calc_median_val(Sign_pairs, "Coefficient")
outNetwork = Puc_compatible_network(pairs_df, microbes_df)

# if this returns TRUE, you did everything correct!
#identical(outNetwork[,c("partner1", "partner2")], Microbe_network_precomputed[,c("partner1", "partner2")])
#write.csv (outNetwork,"microbe-networkFile.csv", quote=FALSE)
identical(rownames(outNetwork[,c("partner1", "partner2")]), rownames(Microbe_network_precomputed[,c("partner1", "partner2")]))
dim(outNetwork[,c("partner1", "partner2")])
dim(Microbe_network_precomputed[,c("partner1", "partner2")])


# plot networks
dfNetwork = outNetwork[,c("partner1", "partner2", "combinedCoefficient")]
print(head(dfNetwork))

#?graph_from_data_frame
# http://kateto.net/networks-r-igraph
g = graph_from_data_frame(dfNetwork,directed = F, vertices = NULL)
print(g, e=TRUE, v=TRUE)
# see the edges and vertices
#E(g) ; V(g) ; edge_attr(g) ; vertex_attr(g)

# plots the networks
plothis = induced.subgraph(g, V(g))
#visualization(plothis,node.size=4,node.label=V(g)$name,node.label.color="blue")
plot(plothis)

cluster1 = Identify_subnetworks(outNetwork)
summary(cluster1)

# if this returns TRUE, you did everything correct!
identical(get.edgelist(cluster1), get.edgelist(microbes_mcode_cluster1_precomputed)) # OR #sum(get.adjacency(cluster1) != get.adjacency(microbes_mcode_cluster1_precomputed)) # should return 0
#write_graph(cluster1, "microbes_mcode_cluster1_edges.txt", "ncol")

plot(cluster1, vertex.label=NA, vertex.size=3)

