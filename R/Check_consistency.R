#' Check consistent direction
#'
#' This function selects elements showing consistent direction across datasets
#'
#' @param df, patternToSearch, threshold
#' @return A matrix with elements showing consistent direction across datasets
#' @export
Check_consistency = function(s_df, pattern, Threshold, numbDatasets=2){
  # select the data of interest
  patternData = s_df[,grep(pattern,colnames(s_df)), drop=FALSE]
  head(patternData)

  rows_passing_consistency = c()
  res_pos = apply(patternData, 1, function(x) sum(x > Threshold))
  rows_passing_consistency = c(rows_passing_consistency, names(res_pos[res_pos==numbDatasets]))

  res_neg = apply(patternData, 1, function(x) sum(x < Threshold))
  rows_passing_consistency = c(rows_passing_consistency, names(res_neg[res_neg==numbDatasets]))

  s_df = s_df[rows_passing_consistency, -1]
  return(s_df)
}
