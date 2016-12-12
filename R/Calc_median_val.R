#' Calculate median coefficient or fold change
#'
#' This function calculate median coefficient or fold change for each element
#'
#' @param df, pattern
#' @return A data frame with statistics for each element
#' @export
Calc_median_val = function(df, pattern){
  Colnames = colnames(df)[grep(pattern,colnames(df))]
  interestedData = df[,Colnames,drop=F]
  interestedData = as.matrix(interestedData)
  interestedData = apply(interestedData,2,function(x){as.numeric(as.vector(x))})
  combined = apply(interestedData,1, function(x){round(median(x, na.rm = TRUE), 3)})
  oldColnames = colnames(df)
  df = cbind(df,combined)
  colnames(df) = c(oldColnames, paste("combined", pattern, sep=''))
  print(head(df))
  return(df)
}
