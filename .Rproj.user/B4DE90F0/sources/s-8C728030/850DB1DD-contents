library(RERconverge)
source("../phenotype_tree_processing.R")


getPathvecfromForegroundTree=function(foreground_tree){
  fg_edge_idx = which(foreground_tree$edge.length==1)
  fg_node_idx = foreground_tree$edge[fg_edge_idx,2]
  foreground_species = foreground_tree$tip.label[fg_node_idx]
  
  pathvec = rep(0, length(foreground_tree$tip.label))
  names(pathvec) = foreground_tree$tip.label
  pathvec[foreground_species] = 1
  pathvec
}


neutraltree = read.tree("../toydatasets/neutraltree_24way.tree")
foregrounds = scan("../toydatasets/subterranean_foregrounds.txt", what="character")
root_species = "mm10"

num_perms = 10


getPermulatedPhenotypeTree=function(foregrounds, neutraltree, root_species){
  foreground_tree = getTreeFromForegrounds(foregrounds, neutraltree, plotTree=F)
  pathvec = getPathvecfromForegroundTree(foreground_tree)
  
  masterTree = list()
  masterTree[[1]] = neutraltree
  names(masterTree) = c("masterTree")
  
  permulated_tree = simBinPhenoSSM(foreground_tree, masterTree, root_species, foregrounds, sisters_list=NULL, pathvec, plotTreeBool = F)
  permulated_tree
}
