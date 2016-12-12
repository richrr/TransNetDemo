#' Get counts for all shortest paths between given nodes in a bipartite networks
#'
#' This function calculates shortest paths between given nodes in a bipartite networks
#'
#' @param pair, myNetwork, sumAllFractionsForAllNodes
#' @return Updates the row in pair,sumAllFractionsForAllNodes
#' @export
Get_shortest_paths = function(pair,myNetwork){  # for each row in "sumAllFractionsForAllNodes", calculate the values for the columns and update this row in "sumAllFractionsForAllNodes"
  #print(pair)
  allShortestPaths = get.all.shortest.paths(myNetwork, from= pair[1], to= pair[2], mode = "all", weights=NULL) # calculate the shortest paths
  #print(allShortestPaths)

  count = data.frame()
  if (length(allShortestPaths$res)==0){  # if there is no shortest path, do not update "sumAllFractionsForAllNodes"
    #fractions = cbind(FromNodes,rep(0,times=length(FromNodes)))
    #print("no paths")
  }else{
    allShortestPaths = do.call(rbind,allShortestPaths$res) # a matrix with each row containing a shortest path between the two nodes
    #print(allShortestPaths)

    nodesInPath = as.vector(allShortestPaths[,c(-1,-ncol(allShortestPaths))]) # get rid of the two nodes that are under study, the nodes in the shortest paths are extracted into this vector
    #print(nodesInPath)
    count = table(as.factor(nodesInPath))/nrow(allShortestPaths) # for each node in the shortest paths, calculate how many times(normalized) it appears
    #print(count)

    v = names(count)
    count = data.frame(rbind((count)))
    colnames(count) = v
    rownames(count) = paste(pair, collapse='_')
    #print(count)
  }
  count
}
