# GeneExpression
## Purpose
Takes a bam file and annotation file, generates a table of FPKM for each gene
## Requirements
Written in R.

Requires:

[argparse](https://cran.r-project.org/web/packages/argparse/index.html)

[Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html)

[DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)

[GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)

[GenomicAlignments](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html)

[BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html)

# Usage
Basic usage is:
```bash
GeneExpression.R -b bamfile1.bam,bamfile2.bam -s sample1,sample2 -a annotation.gff
```
GeneExpression accepts the following additional arguments (all of which have defaults already set):

-c (number of cores)
