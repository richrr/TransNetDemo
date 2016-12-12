#' Get_template_matrix
#'
#' This function creates a matrix , with each entry recording for each node(in a column),
#' if it appears(not 1 or 0, but a normalized number, since there can be more
#' than one shortest path between a pair of nodes) in the shortest path between a
#' pair of nodes(in a row), the rows contain all pairs of nodes between two groups.
#' Thus the sum of each column is the betweeness centrality for the node(in that column)
#'
#' @param FromNodes, allPairs
#' @return A matrix with rownames (source-target) and column names as nodes from source nodes set.
#' @export
Get_template_matrix = function(FromNodes, allPairs){
  sumAllFractionsForAllNodes = matrix(0,ncol=length(FromNodes)
                                      ,nrow=nrow(allPairs)
                                      ,dimnames= list(
                                        paste(as.vector(t(allPairs[,1])),as.vector(t(allPairs[,2])),sep="_")
                                        , as.character(FromNodes)
                                      )
  )
  return(sumAllFractionsForAllNodes)
}
