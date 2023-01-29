library(optparse)

option_list = list(make_option(c('-u', '--fgunmerged'), type='character', default=NULL, help='path to BED file containing the coordinates of foreground elements (unmerged)', metavar='character'),
		   make_option(c('-t', '--motif_file'), type='character', default=NULL, help='path to BED file containing the coordinates of calls for a specific motif', metavar='character'),
		   make_option(c('-o', '--output_folder'), type='character', default=NULL, help='path to output folder for the specific motif', metavar='character'),
		   make_option(c('-x', '--prefix'), type='character', default=NULL, help='character prefix for the output files for this motif', metavar='character')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

fg_file = opt$fgunmerged
motif_file = opt$motif_file
output_folder = opt$output_folder
prefix = opt$prefix

system(paste('mkdir', output_folder))

output_path_motif = paste0(output_folder, prefix, 'motif_coords.bed')
output_path_element = paste0(output_folder, prefix, 'element_coords.bed')

command_text_motif = paste0('bedtools intersect -a ', motif_file, ' -b ', fg_file, ' -u > ', output_path_motif)
system(command_text_motif)

command_text_element = paste0('bedtools intersect -a ', fg_file, ' -b ', motif_file, ' -u > ', output_path_element)
system(command_text_element)

