library(RERconverge)
source("../phenotype_tree_processing.R")
source("../permulation_functions.R")


neutraltree = read.tree("../toydatasets/neutraltree_24way.tree")
foregrounds = scan("../toydatasets/subterranean_foregrounds.txt", what="character")
root_species = "mm10"

num_perms = 10

test = getPermulatedPhenotypes(foregrounds, neutraltree, 20, "mm10")

