library(optparse)
source("functions/loadPhyloConverge.R")
require("doParallel")
require("foreach")

option_list = list(make_option(c("-n", "--neutralmodpath"), type="character", default=NULL, help="file (.mod) path to neutral evolution model", metavar="character"),
		   make_option(c("-p", "--perm"), type="character", default=NULL, help="file (.RDS) path to the list of permulated phenotypes", metavar="character"),
		   make_option(c("-f", "--foregrounds"), type="character", default=NULL, help="txt file containing the names of the foreground species", metavar="character")
		   make_option(c("-a", "--alnmasterfolder"), type="character", default=NULL, help="master folder containing all the alignments to score", metavar="character"),
		   make_option(c("-o", "--outfolder"), type="character", default="outfolder/", help="output folder", metavar="character"),
		   make_option(c("-r", "--refseq"), type="character", default=NULL, help="reference species", metavar="character"),
		   make_option(c("-c", "--numcores"), type="integer", default=2, help="number of cores to parallelize over", metavar="integer"))

print('Parsing arguments')
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

neutralmod_path = opt$neutralmodpath
permulated_foregrounds_path = opt$perm
aln_master_folder = opt$alnmasterfolder
outfolder = opt$outfolder
refseq = opt$refseq
no_cores = opt$numcores
foregrounds_path = opt$foregrounds


print('Loading inputs')
neutralMod = read.tm(neutralmod_path)
permulated_foregrounds = readRDS(permulated_foregrounds_path)
foregrounds = readLines(foregrounds_path)

print('Exporting objects to cores')
clust = makeCluster(no_cores)
clusterExport(clust, "neutralMod")
clusterExport(clust, "foregrounds")
clusterExport(clust, "permulated_foregrounds")
clusterExport(clust, "aln_master_folder")
clusterExport(clust, "outfolder")
clusterExport(clust, "refseq")

registerDoParallel(clust)


all_alns = list.files(aln_master_folder)

phyloConvergeForParallel=function(aln){
	source("functions/loadPhyloConverge.R")
	require("doParallel")
	require("foreach")
		
	maf = read.msa(paste0(aln_master_folder, aln))
	score = phyloConverge(foregrounds, permulated_foregrounds, neutralMod, maf, refseq)
	savePath = paste0(outfolder, sub('.fa', '_scores.RDS', aln))
	saveRDS(score, savePath)
}

print('Run parallel')
parLapply(clust, all_alns, phyloConvergeForParallel)
stopCluster(clust)

