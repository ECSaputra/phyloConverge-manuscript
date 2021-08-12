
#'Converts a binary phenotype tree to a vector of foreground species names
#' @param tree a phylo object with binary edge lengths (foreground edges = 1, background edges = 0)
#' @return a vector of foreground species names
#' @export
getForegroundsFromTree=function(tree){
  edge_length = tree$edge.length
  edge_nodes = tree$edge
  
  foreground_edges = edge_nodes[which(edge_length==1),]
  num_tips = length(tree$tip.label)
  
  tip_foregrounds = tree$tip.label[foreground_edges[which(foreground_edges[,2] <= num_tips),2]]
  ind_internal_foregrounds = which(foreground_edges[,2] > num_tips)
  if (length(ind_internal_foregrounds) > 0){
    internal_foregrounds = NULL
    for (i in 1:length(ind_internal_foregrounds)){
      foreground_branch = foreground_edges[ind_internal_foregrounds[i],]
      internal_foregrounds = c(internal_foregrounds, getInternalBranchName(foreground_branch[2], tree, foreground_edges, num_tips))
    }
    foregrounds = c(tip_foregrounds, internal_foregrounds)
  } else {
    foregrounds = tip_foregrounds
  }
  foregrounds
}

#'Produces phyloConverge-ready name of an ancestral branch
#' @param node internal node ID number, corresponding to the node ID names in the tree object
#' @param tree a phylo object with the tree topology of interest
#' @param allEdges a matrix containing all edges (e.g., the 'edge' component in the tree object, or a subset of it)
#' @param num_tips total number of tip species in the topology
#' @return branch_name the name of the input ancestral node
#' @export
getInternalBranchName=function(node,tree,allEdges,num_tips){
  daughter_nodes = getTipDaughterNodes(node,allEdges,num_tips)
  reprDaughters = tree$tip.label[c(min(daughter_nodes), max(daughter_nodes))]
  branch_name = paste(reprDaughters, collapse="-")
  branch_name
}

#' @keywords internal
getDaughterNodes=function(node, allEdges){
  daughters = allEdges[which(allEdges[,1] == node),2]
  daughters
}

#' @keywords internal
getTipDaughterNodes=function(node, allEdges, num_tips){
  if (node <= num_tips){
    daughters=NULL
  } else {
    alldaughters=allEdges[which(allEdges[,1]==node),2]
    ind_tip_daughters = which(alldaughters<=num_tips)
    if (length(ind_tip_daughters) == 2){
      daughters=alldaughters
    } else {
      ind_anc_daughters = which(alldaughters > num_tips)
      for (i in 1:length(ind_anc_daughters)){
        daughters = c(alldaughters[ind_tip_daughters],getTipDaughters(alldaughters[ind_anc_daughters[i]], allEdges, num_tips))
      }
    }
  }
  daughters
}

#' @keywords internal
getTipDaughterNames=function(node,tree,allEdges,num_tips){
  daughter_nodes = getTipDaughterNodes(node,allEdges,num_tips)
  daughter_names = tree$tip.label[sort(daughter_nodes, decreasing=F)]
  daughter_names
}

