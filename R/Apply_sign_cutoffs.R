#' Select significant elements
#'
#' This function applies the following significance cutoffs
#' individual p-value < 0.3; combined p-value < 0.05; fdr < 0.1
#'
#' @param df, individualPvalueCutoff, combinedPvalueCutoff, combinedFDRCutoff
#' @return A matrix with statistics for signifcant elements
#' @export
Apply_sign_cutoffs = function(df, individualPvalueCutoff = 0.3, combinedPvalueCutoff = 0.05, combinedFDRCutoff = 0.1){
  # find significant in individual pvalue: all of the pvalues for all of the datasets for each element must be smaller than threshold
  # find pvalue data
  # if it is a single dataset it becomes a vector so adding the drop=FALSE
  pvalueData = df[,grep("pvalue",colnames(df)), drop=FALSE]
  print(head(pvalueData))
  pvalueData = as.matrix(pvalueData)
  pvalueData = apply(pvalueData,2,function(x){as.numeric(as.vector(x))})

  # calculate the largest pvalue among all datasets for each gene, this largest pvalue must be smaller than threshold
  passIndevidualPvalue = apply(pvalueData,1,max)<individualPvalueCutoff
  Sign_elements = df[passIndevidualPvalue,]
  Sign_elements <- Sign_elements[Sign_elements$"combinedPvalue" < combinedPvalueCutoff  &  Sign_elements$"combinedFDR" < combinedFDRCutoff, ]
  return(Sign_elements)
}
