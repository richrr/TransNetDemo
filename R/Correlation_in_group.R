#' Correlation_in_group
#'
#' This function calculates correlation in a group for each element pair
#'
#' @param mapFile, geneFile, sampleIdColName, factorColName, groupA, geneSymbolColName, selected genes
#' @return A data frame with statistics for each element
#' @export
Correlation_in_group = function(mapFile, geneFile, sampleIdColName, factorColName, groupA,  geneSymbolColName, selected_genes, usepairs=NA){

  ########################### extract inputs ###############################
  # read the mapping file
  map = read.delim(mapFile, header=T)
  rownames(map) = map[,sampleIdColName]
  head(map)

  # select samples from group
  hfhs_samples = map[which(map[,factorColName] == groupA),]
  head(hfhs_samples)

  # read the data file
  genes = read.delim(geneFile, header=T, check.names = F)
  rownames(genes) = genes[,geneSymbolColName]
  head(genes)

  # create gene pairs if not already provided
  pairs = ''
  if(is.na(usepairs)) {
    pairs = t(combn(selected_genes,2))[,2:1]
  }else{
    pairs = usepairs
  }
  head(pairs)


  ############################# calculate correlation for each pair ##################################
  Corr_pairs = apply(pairs, 1, Calc_cor, genes, rownames(hfhs_samples))


  ############################# format output ##################################
  # add the pair as the column name
  colnames(Corr_pairs) = paste(as.vector(pairs[,1]),as.vector(pairs[,2]),sep="<==>")
  # the pairs become row names; and the corr coeff and pvalue are the column names
  Corr_pairs = t(Corr_pairs)
  # append method used to column labels
  colnames(Corr_pairs) = c(paste(groupA, "Coefficient",sep="_"), "pvalue")


  ############################# calculate FDR  #############################
  FDR = p.adjust(Corr_pairs[,colnames(Corr_pairs)[grep("pvalue",colnames(Corr_pairs))]],method="fdr")
  Corr_pairs = as.data.frame(cbind(Corr_pairs, FDR))

  return(Corr_pairs)
}
