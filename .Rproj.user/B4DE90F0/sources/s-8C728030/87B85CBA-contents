library(RERconverge)
source("phenotype_tree_processing.R")

getPermulatedPhenotypeTree=function(foregrounds, neutraltree, root_species){
  foreground_tree = getTreeFromForegrounds(foregrounds, neutraltree, plotTree=F)
  pathvec = getPathvecfromForegroundTree(foreground_tree)
  
  masterTree = list()
  masterTree[[1]] = neutraltree
  names(masterTree) = c("masterTree")
  
  permulated_tree = simBinPhenoSSM(foreground_tree, masterTree, root_species, foregrounds, sisters_list=NULL, pathvec, plotTreeBool = F)
  permulated_tree
}
