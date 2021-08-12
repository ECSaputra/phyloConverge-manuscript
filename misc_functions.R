#'Subset foreground species based on existing species in alignment file
#' @param align alignment object read from read.msa
#' @param foregrounds character vector containing the foreground species
#' @return foregrounds.out character vector containing the edited foregrounds
#' @export
checkForegrounds=function(align, foregrounds){
	species.in.alignment = align$names
        tip.foregrounds = foregrounds[which(!grepl("-", foregrounds))]
        foregrounds.out = tip.foregrounds[which(tip.foregrounds %in% species.in.alignment)]

	internal.nodes = foregrounds[which(grepl("-", foregrounds))]
	if (length(internal.nodes) > 0){
		for (i in 1:length(internal.nodes)){
			internal.node = internal.nodes[i]
			daughter.tips = strsplit(internal.node,"-")[[1]]
			
			if (length(which(daughter.tips %in% species.in.alignment)) == 2){
				foregrounds.out = c(foregrounds.out, internal.node)
			}
		}
	}
	foregrounds.out
}
