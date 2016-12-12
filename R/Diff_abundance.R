#' Compare abundance between groups for an element
#'
#' This function compares the values in two groups for an element
#'
#' @param element, df, samplesA, samplesB
#' @return A matrix with statistics for each element
#' @export
Diff_abundance = function(element, df, samplesA, samplesB){
  samplesA = as.vector(samplesA)
  samplesB = as.vector(samplesB)

  Vec1 = as.numeric(df[element, samplesA])
  Vec2 = as.numeric(df[element, samplesB])

  meanVec1 = mean(Vec1)
  meanVec2 = mean(Vec2)
  fold_change_mean = meanVec1/meanVec2

  medianVec1 = median(Vec1)
  medianVec2 = median(Vec2)
  fold_change_median = medianVec1/medianVec2

  p = t.test(Vec1, Vec2)

  Result = as.matrix(c(meanVec1, meanVec2, fold_change_mean, medianVec1 , medianVec2 , fold_change_median, p$p.value))

  Result
}
