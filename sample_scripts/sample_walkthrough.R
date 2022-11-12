source('functions/loadPhyloConverge.R')

'''
getting neutral tree, foreground tree, etc
'''

neutralMod_path = '../data/eyeStudy.tree.mod'
foregrounds = c('nanGal1', 'chrAsi1', 'conCri1', 'hetGla2')

neutral_tree = getTreeFromEvolMod(neutralMod_path)

foreground_tree = getTreeFromForegrounds(foregrounds, neutral_tree)
foreground_names = getForegroundsFromTree(foreground_tree)

'''
getting permulated phenotypes
'''

num_perms = 100
root_species = "mm10"

permulated_phenotype_trees = getPermulatedPhenotypes(foregrounds, neutral_tree, num_perms, root_species, output_mod="trees")
permulated_foregrounds = getPermulatedPhenotypes(foregrounds, neutral_tree, num_perms, root_species, output_mod="names")


'''
running phyloConverge
'''

### read the neutral evolution model and the alignment using RPHAST functions
alnfile = '../data/CNEs-alignment/chr1/CNE027703.fa'

msa = read.msa(alnfile)
neutralMod = read.tm(neutralMod_path)

### scoring the entire alignment: leave feature as NULL (default)
score_whole_adaptive = phyloConverge(foregrounds, permulated_foregrounds, neutralMod, msa, "mm10")
score_whole_nonadp = phyloConverge(foregrounds, permulated_foregrounds, neutralMod, msa, "mm10", adapt=F)


### scoring a specific region in the alignment
bed = data.frame("chr"="chr1", "start"=5, "end"=20, "name"="feature1")
feature = convertBedToFeature(bed, "mm10")
score_subregion = phyloConverge(foregrounds, permulated_foregrounds, neutralMod, msa, "mm10", feature, adapt=F)

### scoring multiple subregions in the alignment
bed = data.frame("chr"=rep("chr1",3), "start"=c(5, 25, 50), "end"=c(20, 35, 100), "name"=paste0("feature", c(1,2,3)))
features = convertBedToFeature(bed, "mm10")

score_multiple = phyloConverge_new(foregrounds, permulated_foregrounds, neutralMod, msa, "mm10", features, adapt=T)

print('test')