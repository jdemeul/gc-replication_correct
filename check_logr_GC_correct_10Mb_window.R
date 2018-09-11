### check effect of CNAs on 10Mb correlation window

library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)

CNDIR <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/"

sampleid <- "fc9d93b6-92e8-acb7-e040-11ac0d487dee"

cnfile <- file.path(CNDIR, paste0(sampleid, ".consensus.20170119.somatic.cna.annotated.txt"))
cnas <- read.delim(file = cnfile, as.is = T)

cnas_gr <- GRanges(seqnames = cnas$chromosome, ranges = IRanges(start = cnas$start, end = cnas$end),
                   total_cn = cnas$total_cn, seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))

# debug(get_gcCorrected_logr)

correctiondata <- read_gc_content_replic_timing(gc_content_file_prefix = GCCORRECTPREFIX,
                                                replic_timing_file_prefix = REPCORRECTPREFIX,
                                                chrom_names = chrom_names)

TUMOURNAME <- sampleid
chrom_names <- c(1:22, "X")

Tumour_LogR_file=paste(TUMOURNAME,"_mutantLogR.tab", sep="")
outfile=paste(TUMOURNAME,"_mutantLogR_gcCorrected.tab", sep="")
correlations_outfile=paste(TUMOURNAME, "_GCwindowCorrelations.txt", sep="")
gc_content_file_prefix=NULL
replic_timing_file_prefix=NULL
recalc_corr_afterwards = T
chrom_names=chrom_names
gc_content_df=correctiondata$gc
replic_timing_df=correctiondata$rep


