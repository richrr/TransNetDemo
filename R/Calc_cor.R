#' Calculate correlation between vectors
#'
#' This function calculate correlation between vectors
#'
#' @param pair, df, samples
#' @return A matrix with statistics for each element pair
#' @export
Calc_cor = function(pair, df, samples){
  idxs = as.vector(samples)

  c1 = as.numeric(df[pair[1], idxs])
  c2 = as.numeric(df[pair[2], idxs])

  p = cor.test(c1,c2)
  outLine = as.matrix(c(p$estimate,p$p.value))
  outLine
}

