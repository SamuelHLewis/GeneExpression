#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("BiocParallel"))

#########################
# user argument parsing #
#########################

# create parser object
parser <- ArgumentParser()
# parse bam file names
parser$add_argument("-b", "--bamfiles", type="character", help="Bam files")
# parse sample names
parser$add_argument("-s", "--samplenames", type="character", help="Sample names corresponding to each bam file")
# parse annotation file names
parser$add_argument("-a", "--annotationfile", type="character", help="Annotation file")
# parse cores
parser$add_argument("-c", "--cores", type="character", help="Number of cores to use")
# collect arguments
args <- parser$parse_args()
# assign arguments to variables
baminput = strsplit(args$bamfiles,",")[[1]]
samples = strsplit(args$samplenames,",")[[1]]
annotation = args$annotationfile
# check that the same number of bam files and sample names have been given
if (length(baminput)==length(samples)) {
        sprintf("Each bam file has an accompanying sample name - looking good")
} else {
        stop("ERROR: there are different numbers of bam files and sample names")
}

register(MulticoreParam(workers = args$cores))

# create a bamfile list
bamfiles=BamFileList(baminput)
# check information for first file
seqinfo(bamfiles[1])
# read annotation file and make it into a TxDb object
gfffile=file.path(annotation)
txdb=makeTxDbFromGFF(gfffile,format="gff",circ_seqs=character())
txdb
# make a GRangesList of all exons grouped by gene (each element is a list of exons for a given gene)
ebg=exonsBy(txdb,by="gene")
ebg
# create a summarized experiment object, which generates counts for each gene
se=summarizeOverlaps(features=ebg, reads=bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=FALSE, fragments=TRUE)
# add the sample names to the summarized experiment object
colData(se)<-DataFrame(baminput,samples)
# check that the se object contains necessary information
colData(se)
# create a DESeqDataSet object for the counts for each gene, specifying the design as baminput by sample (NB: this is irrelevant for just plotting FPKM, but a required argument to make the object)
dds=DESeqDataSet(se, design = ~ 1)
# calculate FPKM values for each gene
scaledcounts=fpkm(dds,robust=TRUE)
# output FPKM to file
scaledcountsDF=as.data.frame(scaledcounts)
colnames(scaledcountsDF)=samples
write.csv(scaledcountsDF,file="FPKM.csv")

