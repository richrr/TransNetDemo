#' Calculate network that satisfies the PUC criteria
#'
#' This function keeps edges that satisfy expected correlation and fold change relationships
#'
#'
#' @param pairs_df, genes_df
#' @return A data frame (network) that is Puc compatiable
#' @export
Puc_compatible_network = function(pairs_df, genes_df){
  #--------------------------------------------------------------------------------
  # calculate PUC
  #--------------------------------------------------------------------------------

  #change out format to partner1 partner2 for PUC
  row_names_pairs_df = rownames(pairs_df)
  head(row_names_pairs_df)
  pair = stringr::str_split( row_names_pairs_df ,"<==>")
  pairs = t(as.data.frame(pair))

  colnames(pairs) = c("partner1","partner2")
  pairs_df = apply(pairs_df, 2, function(x) as.numeric(as.character(x))) # convert the chars to numeric
  rownames(pairs) = row_names_pairs_df # remove this if you do not want row names
  outForPUC = cbind(pairs,pairs_df)
  head(outForPUC)


  # select pvalue columns
  grep_cols_c = grep("pvalue", colnames(outForPUC), value=TRUE, fixed=TRUE)
  # select "combinedPvalue", "combinedFDR", "combinedCoefficient"
  grep_cols_c = append(grep_cols_c, grep("combined" , colnames(outForPUC), value=TRUE, fixed=TRUE) )
  # select the partner columns
  g_grep_cols = c("partner1","partner2", grep_cols_c)

  outForPUC = outForPUC[,g_grep_cols]
  head(outForPUC)


  # attach the foldChange information for each partner
  FoldChangeCol = grep(paste("combined" , foldchVar, sep=''), colnames(genes_df), value=TRUE, fixed=TRUE)
  FoldMetab1_InPair = genes_df[as.vector(outForPUC[,"partner1"]), FoldChangeCol, drop=F]
  colnames(FoldMetab1_InPair) = c("partner1_FoldChange")
  FoldMetab2_InPair = genes_df[as.vector(outForPUC[,"partner2"]), FoldChangeCol, drop=F]
  colnames(FoldMetab2_InPair) = c("partner2_FoldChange")
  outForPUC = cbind(outForPUC,FoldMetab1_InPair,FoldMetab2_InPair)
  head(outForPUC)

  # calculate correlation Direction For combined correlation coefficient of interest
  # at this point we only have the consistent pairs left, so the value of combined corr coeff is ok to use
  interestedCoefficientColnames = grep("Coefficient",colnames(outForPUC), value=TRUE, fixed=TRUE)
  print(interestedCoefficientColnames)
  interestedCorrelationData = outForPUC[,interestedCoefficientColnames, drop=FALSE]
  interestedCorrelationData = apply(interestedCorrelationData,2,function(x){as.numeric(as.vector(x))})
  correlationDirection = interestedCorrelationData/abs(interestedCorrelationData)


  # calculate fold change direction for each partner
  FoldChangeColnames = colnames(outForPUC)[grep("FoldChange",colnames(outForPUC))] # since this is using the combined fold change calculated above, you do not need the foldchVar variable

  FoldChangeData = outForPUC[,FoldChangeColnames]
  FoldChangeDirection = (FoldChangeData-1)/abs(FoldChangeData-1)
  names(FoldChangeDirection) = c()
  colnames(FoldChangeDirection) = paste(colnames(FoldChangeData),"Direction",sep="_")

  # calculate if fold change direction are the same for the two partners
  IfFoldChangeDirectionMatch = apply(FoldChangeDirection,1,prod)
  names(IfFoldChangeDirectionMatch) = c()
  colnames(correlationDirection) = c()

  # use "correlationDirection" and "IfFoldChangeDirectionMatch" to calc PUC,
  # i.e. if these two are the same PUC=1 (good)
  PUC = IfFoldChangeDirectionMatch * correlationDirection
  outForPUC = cbind(outForPUC,correlationDirection,FoldChangeDirection,IfFoldChangeDirectionMatch,PUC)
  head(outForPUC)

  #write.csv(outForPUC, "PUC-output.csv" ,row.names=FALSE)

  # find PUC expected
  out = outForPUC[outForPUC[,"PUC"]==1 ,]
  outNetwork = out[!is.na(out$"PUC"),] # remove the rows with 'NA' in PUC columns
  return(outNetwork)

}
