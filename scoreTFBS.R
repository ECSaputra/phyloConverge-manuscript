library(rphast)
library(phyloConverge)
library(RERconverge)
library(optparse)

option_list = list(make_option(c('-p', '--permulated_foregrounds_path'), type='character', default=NULL, help='RDS file path to a list object containing permulated phenotypes', metavar='character'),
		   make_option(c('-t', '--treepath'), type='character', default=NULL, help='path to mod file representing the neutral substitution model', metavar='character'),
		   make_option(c('-f', '--foregrounds'), type='character', default=NULL, help='foregrounds species names, split by commas', metavar='character'),
		   make_option(c('-w', '--aln_folder'), type='character', default=NULL, help='folder path containing alignments for each element', metavar='character'),
		   make_option(c('-r', '--refseq'), type='character', default=NULL, help='reference species', metavar='character'),
		   make_option(c('-i', '--int_folder'), type='character', default=NULL, help='folder path to intersecting motifs info', metavar='character'),
		   make_option(c('-o', '--output_path'), type='character', default=NULL, help='RDS output file path', metavar='character')
)



opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

permulated_foregrounds_path = opt$permulated_foregrounds_path
treepath = opt$treepath
foregrounds_text = opt$foregrounds
aln_folder = opt$aln_folder
refseq = opt$refseq
int_folder = opt$int_folder
output_path = opt$output_path


permulated_foregrounds = readRDS(permulated_foregrounds_path)
neutralMod = read.tm(treepath)
foregrounds = strsplit(foregrounds_text, ',')[[1]]

int_folder_content = list.files(int_folder)

motif_coords_file = int_folder_content[which(grepl('motif_coords', int_folder_content))]
element_coords_file = int_folder_content[which(grepl('element_coords', int_folder_content))]

motif_coords = read.delim(paste0(int_folder, motif_coords_file), header=F)
elements_coords = read.delim(paste0(int_folder, element_coords_file), header=F)

motif_counter = 0
motif_scores = NULL
for (ii in 1:nrow(elements_coords)){
	aln = read.msa(paste0(aln_folder, elements_coords[ii,4], '.fa'))
	motifs_in_chr = motif_coords[which(motif_coords[,1]==elements_coords[ii,1]),]
	if (aln$offset==0){
		aln$offset = elements_coords[ii,2]-1
	}
	motifs_in_aln = getElementsInMaf(motifs_in_chr, aln)
	if (!is.null(motifs_in_aln)){
		motif_counter = motif_counter + nrow(motifs_in_aln)
		features = convertBedToFeature(motifs_in_aln, refseq)
		motif_scores_in_aln = phyloConverge(foregrounds, permulated_foregrounds, neutralMod, aln, refseq, features, adapt=T)
		colnames(motifs_in_aln) = c('chr', 'start', 'end', 'name')
		motifs_in_aln$permPval = motif_scores_in_aln$permPval
		motifs_in_aln$corr_score = motif_scores_in_aln$corr_score
		motifs_in_aln$uncorr_score = motif_scores_in_aln$uncorr_score
		motif_scores = rbind(motif_scores, motifs_in_aln)
	}
}

saveRDS(motif_scores, output_path)

