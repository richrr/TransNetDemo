#library(TransNetDemo)
library(gplots)


########################### extract inputs - study 1 ###############################

# groups to compare
groupA = "HFHS"
groupB = "NCD"
sampleIdColName = "SampleID"
factorColName = "Factor"
geneSymbolColName ="IdSymbol"

# map files
mapf1 = system.file("extdata", "mapping_file.rand.1.tsv", package = "TransNetDemo")

# read the mapping file
map1 = read.delim(mapf1, header=T)
rownames(map1) = map1[,sampleIdColName]

# select samples from groups
hfhs_samples1 = as.vector(map1[which(map1[,factorColName] == groupA),sampleIdColName])
ncd_samples1 = as.vector(map1[which(map1[,factorColName] == groupB),sampleIdColName])

# colors as per the group
samplecolors <- c(rep("magenta",25) , rep("blue",25))

# gene files
genef1 = system.file("extdata", "gene_file_1.tsv", package = "TransNetDemo")

# read the data file
genes1 = read.delim(genef1, header=T, check.names = F, row.names=1)

# make the heat map
dataset1=genes1[rownames(Sign_genes_precomputed_1), c(hfhs_samples1, ncd_samples1)]
dim(dataset1)
heatmap.2(as.matrix(dataset1), col=redgreen(75), ColSideColors=samplecolors, scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)


# microbe files
microbef1 = system.file("extdata", "microbe_file_1.tsv", package = "TransNetDemo")

# read the data file
microbes1 = read.delim(microbef1, header=T, check.names = F, row.names=1)

# make the heat map
dataset1=microbes1[rownames(Sign_microbes_precomputed_1), c(hfhs_samples1, ncd_samples1)]
dim(dataset1)
heatmap.2(as.matrix(dataset1), col=redgreen(75), ColSideColors=samplecolors, scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
