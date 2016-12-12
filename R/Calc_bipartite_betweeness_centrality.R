#' Calc_bipartite_betweeness_centrality
#'
#' This function calculates bipartite_betweeness_centrality for each given element pair
#'
#' @param allPairs, FromNodes, myNetwork
#' @return A data frame with statistics for each element
#' @export
Calc_bipartite_betweeness_centrality = function(allPairs, FromNodes, myNetwork){

  sumAllFractionsForAllNodes = Get_template_matrix(FromNodes, allPairs)
  sumAllFractionsForAllNodes = as.data.frame(sumAllFractionsForAllNodes)
  #print(sumAllFractionsForAllNodes)
  counts = apply(as.matrix(allPairs),1,Get_shortest_paths,myNetwork)
  #print(counts)


   z = lapply(counts, function(x){
     if(nrow(x)>0){
       #print(x)
       columnsForUpdate = as.character(intersect(colnames(x),FromNodes))
       #print(columnsForUpdate)
       rowsForUpdate = rownames(x)
       #print(rowsForUpdate)
       #print(x[1, columnsForUpdate])
       sumAllFractionsForAllNodes[rowsForUpdate,columnsForUpdate] <<- x[1,columnsForUpdate]
     }
   })

  return(sumAllFractionsForAllNodes)
}

