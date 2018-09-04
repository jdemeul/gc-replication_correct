## fix "broken" GC correction
library(readr)

ids <- c("28e81540-4744-4865-b627-c7c9d8a3c2b8", "d182b67c-c622-11e3-bf01-24c6515278c0", "fc9d93b6-92e8-acb7-e040-11ac0d487dee")

BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2018_GC_correction/data/battenberg_logr"
GCCORRECTPREFIX <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_gcContent/1000_genomes_GC_corr_chr_"
REPCORRECTPREFIX <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_repliTiming/1000_genomes_replication_timing_chr_"
chrom_names <- c(1:22, "X")

source("/srv/shared/vanloo/home/jdemeul/projects/2018_GC_correction/gc-replication_correct/prepare_wgs.R")

### read GC content and replicaiton timing data just once

read_gc_content_replic_timing <- function(gc_content_file_prefix, replic_timing_file_prefix, chrom_names) {
  chrom_idx <- 1:length(chrom_names)
  gc_files <- paste0(gc_content_file_prefix, chrom_idx, ".txt.gz")
  gc_content_df <- do.call(rbind, lapply(gc_files, readr::read_tsv, skip = 1, col_names = F, col_types = "-cinnnnnnnnnnnn------"))
  colnames(gc_content_df) <- c("chr", "Position", paste0(c(25,50,100,200,500), "bp"),
                               paste0(c(1,2,5,10,20,50,100), "kb"))
  
  replic_files <- paste0(replic_timing_file_prefix, chrom_idx, ".txt.gz")
  replic_timing_df <- do.call(rbind, lapply(replic_files, readr::read_tsv, col_types = "cin"))
  return(list(gc = gc_content_df, rep = replic_timing_df))
}

correctiondata <- read_gc_content_replic_timing(gc_content_file_prefix = GCCORRECTPREFIX,
                                                replic_timing_file_prefix = REPCORRECTPREFIX,
                                                chrom_names = chrom_names)

for (id in ids) {
TUMOURNAME <- id

RUN_DIR <- file.path(BASEOUT, paste0(TUMOURNAME, "_allelecounts"))
setwd(RUN_DIR)

# debug(gc.replic.correct.wgs)
gc.replic.correct.wgs(Tumour_LogR_file=paste(TUMOURNAME,"_mutantLogR.tab", sep=""),
                      outfile=paste(TUMOURNAME,"_mutantLogR_gcCorrected.tab", sep=""),
                      correlations_outfile=paste(TUMOURNAME, "_GCwindowCorrelations.txt", sep=""),
                      gc_content_file_prefix=NULL,
                      replic_timing_file_prefix=NULL,
                      recalc_corr_afterwards = T,
                      chrom_names=chrom_names, 
                      gc_content_df=correctiondata$gc,
                      replic_timing_df=correctiondata$rep)

}
