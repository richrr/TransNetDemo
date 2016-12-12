#' Calc combined p value across expts.
#'
#' This function code for calc combined p value across datasets.
#'
#' @param data frame
#' @return A matrix with combined pvalue and fdr
#' @export
Calc_combined = function(s_df){
  pvalueData = s_df[,grep("pvalue",colnames(s_df)), drop=FALSE]
  head(pvalueData)

  total_numb_input_files = ncol(pvalueData)
  interestedPvalueData = as.matrix(pvalueData)
  interestedPvalueData = apply(interestedPvalueData,2,function(x){as.numeric(as.vector(x))})
  combinedPvalue = apply(interestedPvalueData,1
                         ,function(pvalues){
                           pvalues = pvalues[!is.na(pvalues)]
                           statistics = -2*log(prod(pvalues))
                           degreeOfFreedom = 2*length(pvalues)
                           combined = 1-pchisq(statistics,degreeOfFreedom)
                         }
  )
  #calculate FDR for combined pvalue
  combinedFDR = p.adjust(combinedPvalue,method="fdr")
  comb_in_df = cbind(s_df, combinedPvalue, combinedFDR)

  return(comb_in_df)
}
