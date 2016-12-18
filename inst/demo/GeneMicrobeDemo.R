#library(TransNetDemo)
library(stringr)
library(ProNet)
library(igraph)
library(ggplot2)


individualPvalueCutoff = 0.3
combinedPvalueCutoff = 0.05
combinedFDRCutoff = 0.2

# groups to compare
groupA = "HFHS"
groupB = "NCD"
sampleIdColName = "SampleID"
factorColName = "Factor"
gene_microbeSymbolColName ="IdSymbol"
# fold change column based on mean or median
foldchVar = "FoldChange_HFHS_NCD"   # or "FoldChange_Median_HFHS_NCD"

# map files
mapf1 = system.file("extdata", "mapping_file.rand.1.tsv", package = "TransNetDemo")
mapf2 = system.file("extdata", "mapping_file.rand.2.tsv", package = "TransNetDemo")

# gene_microbe files
gene_microbef1 = system.file("extdata", "gene_microbe_file_1.tsv", package = "TransNetDemo")
gene_microbef2 = system.file("extdata", "gene_microbe_file_2.tsv", package = "TransNetDemo")

pairs = expand.grid(V(genes_mcode_cluster1_precomputed)$name,V(microbes_mcode_cluster1_precomputed)$name)

# gene_microbe pairs
Corr_pairs = Correlation_in_group(mapf1, gene_microbef1, sampleIdColName, factorColName, groupA,  gene_microbeSymbolColName, NA, pairs)

# Select significant correlations:
#########  (i) using a significance threshold #########
Sign_pairs = Corr_pairs[which(Corr_pairs$pvalue < 0.05 & Corr_pairs$FDR < 0.1),]
######### ######### ######### #########
# OR
######### (ii) in case you have multiple datasets, we recommend meta-analysis ###########
## we have already run these steps and stored the data under the variable "Sign_gene_microbe_pairs_metaanalysis_precomputed" ##
Corr_pairs1 = Corr_pairs
Corr_pairs2 = Correlation_in_group(mapf2, gene_microbef2, sampleIdColName, factorColName, groupA,  gene_microbeSymbolColName, NA, pairs)
numbDatasets=2  # number of datasets

s_df = merge(Corr_pairs1, Corr_pairs2, by="row.names")
rownames(s_df) = rownames(Corr_pairs1)
s_df = Check_consistency(s_df, "Coefficient", 0, numbDatasets)
comb_in_df = Calc_combined(s_df)
Sign_pairs = Apply_sign_cutoffs(comb_in_df,individualPvalueCutoff, combinedPvalueCutoff, combinedFDRCutoff )

# if this returns TRUE, you did everything correct!
identical(Sign_pairs, Sign_gene_microbe_pairs_metaanalysis_precomputed)
#write.csv (Sign_pairs,"Sign_gene_microbes_pairs_File.csv", quote=FALSE)
#write.csv (Sign_pairs,"gene_microbe-networkFile.csv", quote=FALSE)
######### ######### ######### #########


pairs_df = Calc_median_val(Sign_pairs, "Coefficient")

#change out format to partner1 partner2 for PUC
row_names_pairs_df = rownames(pairs_df)
head(row_names_pairs_df)
pair = stringr::str_split( row_names_pairs_df ,"<==>")
pairs = t(as.data.frame(pair))

colnames(pairs) = c("partner1","partner2")
rownames(pairs) = row_names_pairs_df
outForPUC = cbind(pairs,pairs_df)
head(outForPUC)

dfNetwork = outForPUC[,c("partner1", "partner2")]
head(dfNetwork)

g = graph_from_data_frame(dfNetwork,directed = F, vertices = NULL)

# if this returns TRUE, you did everything correct!
identical(get.edgelist(g), get.edgelist(Gene_Microbe_network_precomputed))
#write_graph(g, "gene_microbe_edges.txt", "ncol")


# create TK (bipartite) network using:
# (i) genes_mcode_cluster1_precomputed
# (ii) microbes_mcode_cluster1_precomputed
# (iii) Gene_Microbe_network_precomputed


net1 = get.edgelist(genes_mcode_cluster1_precomputed)
head(net1)

net2 = get.edgelist(microbes_mcode_cluster1_precomputed)
head(net2)

net3 = get.edgelist(Gene_Microbe_network_precomputed)
head(net3)

TK_Network = graph_from_data_frame(rbind(net1, net2, net3), directed = F)
print(TK_Network, e=TRUE, v=TRUE)

# if this returns TRUE, you did everything correct!
identical(get.edgelist(TK_Network), get.edgelist(TK_Network_precomputed))
#write_graph(TK_Network, "Trans_Kingdom_NetworkFile.txt", "ncol")
#write_graph(TK_Network, "Trans_Kingdom_NetworkFile_indices.txt")


# write the mapping of nodes name to indices and group
nodes = data.frame()
for (vertex in V(TK_Network)) {
  Name = V(TK_Network)$name[vertex]
  Id = vertex-1
  Group = ""
  if(Name %in% as.vector(net1)){
    Group = "gene"
  } else {
    Group = "microbe"
  }
  nodes = rbind(nodes, cbind(Name, Id, Group))
}
colnames(nodes)
#write.table(nodes, "Trans_Kingdom_NetworkFile_nodes.txt", quote=F, row.names = F, sep=' ', col.names = T)

# calc bipartite betweenness centrality
FromNodes = as.numeric(nodes[nodes[,3]=="gene",2])
ToNodes = as.numeric(nodes[nodes[,3]=="microbe",2])
allPairs = expand.grid(FromNodes,ToNodes)
myNetwork = TK_Network
sumAllFractionsForAllNodes = Calc_bipartite_betweeness_centrality(allPairs, FromNodes, myNetwork)

sumAllFractionsForAllNodes

forPlot = colSums(sumAllFractionsForAllNodes)  # calculate betweeness centrality for the nodes in the columns
qplot(forPlot,data=as.data.frame(forPlot) ,geom="density")


TK_Network = set_vertex_attr(TK_Network, "type", index = nodes$Id, as.factor(nodes$Group))
plot(TK_Network, vertex.label = NA, layout=layout_as_tree, vertex.color=c( "pink", "skyblue")[1+(V(TK_Network)$type==1)])


# plot in two rows
#ltypes = c(TRUE, FALSE)[1+(V(TK_Network)$type==1)]
#plot(TK_Network, vertex.label = NA, vertex.size=3, layout=layout_as_bipartite(TK_Network,ltypes, hgap = 500), vertex.color=c( "pink", "skyblue")[1+(V(TK_Network)$type==1)])

