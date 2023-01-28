# phyloConverge analysis of subterranean mammal adaptation

This repository stores supplementary scripts used in "Prediction of local convergent shifts in evolutionary rates with _phyloConverge_ characterizes the phenotypic associations and modularity of regulatory elements" by Saputra et al. (preprint on [_bioRxiv_](https://www.biorxiv.org/content/10.1101/2022.05.02.490345v1)).


## Computing correlations between top-ranking regions with functional data

Correlations between top-ranking regions ('foreground elements') with functional genomics data were computed using the R script `correlateWithFunctionalDatasets.R`. The [BEDtools](https://bedtools.readthedocs.io/en/latest/index.html) suite is required to run this function, and make sure that the binaries are added to the PATH environment variable.

The function takes the following parameters:

* -u, --fgunmerged: path to BED file containing the coordinates of foreground elements (unmerged)

* -m, --fgmerged: path to BED file containing the merged coordinates of foreground elements (obtained by running bedtools merge with using-defined distance parameter -d)

* -v, --validation: path to BED file containing coordinates for the functional genomic dataset used for validation

* -a, -allregionsPath: path to BED file containing the coordinates of all scored elements

* -n, --numperm: number of permutations

* -o, outputpath: RDS output file path
