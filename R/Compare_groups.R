#' Compare groups
#'
#' This function compares the values in two groups for each element
#'
#' @param mapFile, geneFile, sampleIdColName, factorColName, groupA, groupB, geneSymbolColName
#' @return A data frame with statistics for each element
#' @export
Compare_groups = function(mapFile, geneFile, sampleIdColName, factorColName, groupA, groupB, geneSymbolColName){

  ########################### extract inputs ###############################
  # read the mapping file
  map = read.delim(mapFile, header=T)
  rownames(map) = map[,sampleIdColName]
  head(map)

  # select samples from groups
  hfhs_samples = map[which(map[,factorColName] == groupA),]
  head(hfhs_samples)
  ncd_samples = map[which(map[,factorColName] == groupB),]
  head(ncd_samples)

  # read the data file
  genes = read.delim(geneFile, header=T, check.names = F)
  rownames(genes) = genes[,geneSymbolColName]
  head(genes)


  ############################# calc diff abundance ##################################
  out = sapply(rownames(genes), Diff_abundance, genes, rownames(hfhs_samples), rownames(ncd_samples))


  ############################# format output ##################################
  Comp_genes = t(out)
  colnames(Comp_genes) = c(paste("Mean", groupA, sep="_"), paste("Mean", groupB, sep="_"), paste("FoldChange", groupA, groupB, sep="_"),
                     paste("Median", groupA, sep="_"), paste("Median", groupB, sep="_"), paste("FoldChange","Median", groupA, groupB, sep="_"), "pvalue")


  ############################# calculate FDR  #############################
  FDR = p.adjust(Comp_genes[,"pvalue"],method="fdr")
  oldColnames = colnames(Comp_genes)
  Comp_genes = as.data.frame(cbind(Comp_genes, FDR))

  return(Comp_genes)
}
