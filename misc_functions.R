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


#'Identify which genetic elements in a table exist in a given multiple alignment file
#' @param elements_bed a BED format matrix/data frame containing a list of genetic elements with their chromosome and coordinates
#' @param maf an MSA object containing the sequence alignment
#' @return elements_in_maf a BED format matrix/data frame listing the elements that exist in maf
#' @export
getElementsInMaf=function(elements_bed, maf){
	coord_range = coord.range.msa(maf)

	ind_elements_in_maf = intersect(which(elements_bed[,2] >= coord_range[1]), which(elements_bed[,3] <= coord_range[2]))

	if (length(ind_elements_in_maf) > 0){
		elements_in_maf = elements_bed[ind_elements_in_maf,]
	} else {
		elements_in_maf = NULL
	}
	elements_in_maf
}