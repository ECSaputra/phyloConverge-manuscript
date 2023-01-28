library(optparse)

option_list = list(make_option(c("-u", "--fgunmerged"), type="character", default=NULL, help="file path to unmerged bed file for foreground elements", metavar="character"),
		   make_option(c("-m", "--fgmerged"), type="character", default=NULL, help="file path to merged bed file for foreground elements", metavar="character"),
		   make_option(c("-v", "--validation"), type="character", default=NULL, help="file path to bed file for validation dataset", metavar="character"),
		   make_option(c("-a", "--allregionsPath"), type="character", default=NULL, help="file path to bed file for all elements", metavar="character"),
		   make_option(c("-n", "--numperm"), type="integer", default=NULL, help="number of permutations", metavar="integer"),
		   make_option(c("-o", "--outputpath"), type="character", default=NULL, help=".RDS file path to output", metavar="character")
		   )


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

fg_regions_unmerged = opt$fgunmerged
fg_regions_merged = opt$fgmerged
validation_data = opt$validation
allregionsPath = opt$allregionsPath
numperm = opt$numperm
savePath = opt$outputpath


allregions = read.delim(allregionsPath, header=F)
unmerged_fgs = read.delim(fg_regions_unmerged, header=F)

total_num_regions = nrow(allregions)
num_unmerged_fgs = nrow(unmerged_fgs)

# Find foreground regions that intersect the validation dataset regions
system(paste("bedtools intersect -a ", fg_regions_merged," -b ", validation_data," -u > tmp_obs_intersect.bed", sep=""))

info = file.info("tmp_obs_intersect.bed")
if (info$size > 0){
	obs_intersect_regions = read.delim("tmp_obs_intersect.bed", header=F)
	obs_intersect = nrow(obs_intersect_regions)
} else {
	obs_intersect = 0
}

### get nulls
null_intersects = NULL
for (i in 1:numperm){
	print(paste("Running perm",i, "/", numperm))

	### randomly select a matched-sized null set of CNEs
	ind.null = sample(total_num_regions, num_unmerged_fgs, replace=F)
	null.regions = allregions[ind.null,]

	### write into bed files, sort, merge
	write.table(null.regions, "null.regions.bed", quote=F, col.names=F, row.names=F, sep="\t")
	system("bedtools sort -i null.regions.bed > null.regions.sorted.bed")
	system("bedtools merge -i null.regions.sorted.bed -d 50 > null.regions.merged.bed")

	### find intersections
	system(paste("bedtools intersect -a null.regions.merged.bed -b ", validation_data," -u > tmp_null_intersect.bed", sep=""))

	info = file.info("tmp_null_intersect.bed")
	if (info$size > 0){
		null_intersect_regions = read.delim("tmp_null_intersect.bed", header=F)
		null_intersects = c(null_intersects, nrow(null_intersect_regions))
	} else {
		null_intersects = c(null_intersects, 0)
	}
}

Z = (obs_intersect - mean(null_intersects))/sd(null_intersects)

output = list("num_overlaps_obs"=obs_intersect, "num_overlaps_null"=null_intersects, "Z"=Z)

saveRDS(output, savePath)


system("rm tmp_obs_intersect.bed")
system("rm null.regions.bed")
system("rm null.regions.sorted.bed")
system("rm null.regions.merged.bed")
system("rm tmp_null_intersect.bed")


