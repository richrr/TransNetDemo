#' Identify top subnetwork
#'
#' This function returns the induced subgraph of the nodes in the top cluster using mcode
#'
#'
#' @param network in which to identify subnetworks
#' @return A data frame (network) that is Puc compatiable
#' @export
Identify_subnetworks = function(outNetwork){

  # identify subnetworks
  dfNetwork = outNetwork[,c("partner1", "partner2", "combinedCoefficient")]
  print(head(dfNetwork))

  #?graph_from_data_frame
  # http://kateto.net/networks-r-igraph
  g = graph_from_data_frame(dfNetwork,directed = F, vertices = NULL)
  print(g, e=TRUE, v=TRUE)
  # see the edges and vertices
  #E(g) ; V(g) ; edge_attr(g) ; vertex_attr(g)

  # plots the networks
  #cluster(g,method="MCODE",layout="fruchterman.reingold")
  #cluster(g,method="FN",layout="fruchterman.reingold")

  #https://cran.r-project.org/web/packages/ProNet/vignettes/Tutorial.pdf
  # identify clusters
  result <- mcode(g,vwp=0.5,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE)
  summary(result$COMPLEX)

  # plot the top cluster
  cluster1<-induced.subgraph(g,result$COMPLEX[[1]])
  #visualization(cluster1,node.size=4,node.label=V(cluster1)$name,node.label.color="blue")

  return(cluster1)
}
