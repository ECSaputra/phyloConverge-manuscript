# phyloConverge analysis of subterranean mammal adaptation

This repository stores supplementary scripts used in "Prediction of local convergent shifts in evolutionary rates with _phyloConverge_ characterizes the phenotypic associations and modularity of regulatory elements" by Saputra et al. (preprint on [_bioRxiv_](https://www.biorxiv.org/content/10.1101/2022.05.02.490345v1)).

The documentation below describes how to run the functions used to perform the downstream analysis shown in the manuscript. The folder `exdata` contains the example datasets used, and the folder `exoutput` contains the example outputs. 

The phyloConverge software can be installed from [this other repository](https://github.com/ECSaputra/phyloConverge).

## Computing correlations between top-ranking regions with functional data

Correlations between top-ranking regions ('foreground elements') with functional genomics data were computed using the R script `correlateWithFunctionalDatasets.R`. The [BEDtools](https://bedtools.readthedocs.io/en/latest/index.html) suite is required to run this function, and make sure that the binaries are added to the PATH environment variable.

The function takes the following parameters:

* -u, --fgunmerged: path to BED file containing the coordinates of foreground elements (unmerged)

* -m, --fgmerged: path to BED file containing the merged coordinates of foreground elements (obtained by running bedtools merge with using-defined distance parameter -d)

* -v, --validation: path to BED file containing coordinates for the functional genomic dataset used for validation

* -a, --allregionsPath: path to BED file containing the coordinates of all scored elements

* -n, --numperm: number of permutations

* -o, --outputpath: RDS output file path

Below is an example on how to run the function:
```
Rscript correlateWithFunctionalDatasets.R -u exdata/cne-coords/acc_CNEs_phyloConverge.bed -m exdata/cne-coords/acc_CNEs_phyloConverge_merged.bed -v exdata/validation-dataset/eye.specific.mouse.ATAC.E11.5.bed -a exdata/cne-coords/mouseCNEs.bed -n 1000 -o exoutput/correlation-functional-data/correlation_embryonic_eye_specific.RDS
```

The output of the function is a list object containing the following variables:

* num_overlaps_obs: the number of overlaps between the merged foreground elements and the validation dataset

* num_overlaps_null: If numperms = X, a X-length vector containing the numbers of overlaps between merged coordinates of each (null) set of randomly selected elements and the validation dataset

* Z: the Z-score computed from num_overlaps_obs with respect to the distribution of num_overlaps_null


## Computing transcription factor binding site (TFBS)-scale convergence signals

The computation of TFBS-scale convergence scores were performed in two steps. First, the `getIntersectingMotifCoords.R` script was used to get coordinates of motifs that intersect the foreground elements. This script is a wrapper for the bedtools 'intersect' function that is used to output two BED files, containing the coordinates of a specific motif that intersect a set of elements, as well as the coordinates of the intersecting elements, respectively. The script takes the following parameters:

* -u, --fgunmerged: path to BED file containing the coordinates of foreground elements (unmerged)

* -t, --motif_file: path to BED file containing the coordinates of calls for a specific motif

* -o, --output_folder: path to output folder for the specific motif

* -x, --prefix: character prefix for the output files for this motif

Below is an example on how to run the function:
```
Rscript getIntersectingMotifCoords.R -u exdata/cne-coords/acc_CNEs_phyloConverge.bed -t exdata/TFBS-calls-mm10/ZSCA4_HUMAN.H11MO.0.D.mm10.bed -o exoutput/TFBS-analysis/intersecting-motifs/intersect_ZSCA4_HUMAN.H11MO.0.D.mm10/ -x ZSCA4_HUMAN.H11MO.0.D.mm10_
```

In the second step, the `scoreTFBS.R` script was used to compute convergence scores for the identified motif coordinates using phyloConverge. The function takes the following parameters:

* -p, --permulated_foregrounds_path: RDS file path to a list object containing permulated phenotypes

* -t, --treepath: path to mod file representing the neutral substitution model

* -f, --foregrounds: foregrounds species names, split by comma

* -w, --aln_folder: folder path containing alignments for each element

* -r, --refseq: reference species

* -i, --int_folder: folder path to intersecting motifs info

* -o, --output_path: RDS output file path

Below is an example on how to run the function:
```
Rscript scoreTFBS.R -p exdata/permulated-foregrounds/permulated_phenotypes_500perms.RDS -t exdata/tree-model/eyeStudy.tree.mod -f chrAsi1,nanGal1,conCri1,hetGla2 -w exdata/cne-alignments/ -r mm10 -i exoutput/TFBS-analysis/intersecting-motifs/intersect_ZSCA4_HUMAN.H11MO.0.D.mm10/ -o exoutput/TFBS-analysis/TFBS-scores/scores_ZSCA4_HUMAN.H11MO.0.D.mm10.RDS
```

