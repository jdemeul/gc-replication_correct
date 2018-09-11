## use Battenberg to get GC-corrected LogR data
library(readr)
library(rslurm)

RELEASETABLEFILE <- "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv"
GCCORRECTPREFIX <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_gcContent/1000_genomes_GC_corr_chr_"
REPCORRECTPREFIX <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_repliTiming/1000_genomes_replication_timing_chr_"

BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/data/battenberg_logR/"

chrom_names <- c(1:22, "X")

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/genecon/GC_utils.R")
source("/srv/shared/vanloo/home/jdemeul/projects/2018_GC_correction/gc-replication_correct/prepare_wgs.R")

### load one time only
releasetable <- read_pcawg_release_table(release_table_file = RELEASETABLEFILE)

### read GC content and replication timing data just once
read_gc_content_replic_timing <- function(gc_content_file_prefix, replic_timing_file_prefix, chrom_names) {
  chrom_idx <- 1:length(chrom_names)
  gc_files <- paste0(gc_content_file_prefix, chrom_idx, ".txt.gz")
  gc_content_df <- do.call(rbind, lapply(gc_files, readr::read_tsv, skip = 1, col_names = F, col_types = "-cinnnnnnnnnnnnnnnnnn"))
  colnames(gc_content_df) <- c("chr", "Position", paste0(c(25,50,100,200,500), "bp"),
                               paste0(c(1,2,5,10,20,50,100,200,500), "kb"),
                               paste0(c(1,2,5,10), "Mb"))
  
  replic_files <- paste0(replic_timing_file_prefix, chrom_idx, ".txt.gz")
  replic_timing_df <- do.call(rbind, lapply(replic_files, readr::read_tsv, col_types = "cin"))
  return(list(gc = gc_content_df, rep = replic_timing_df))
}

correctiondata <- read_gc_content_replic_timing(gc_content_file_prefix = GCCORRECTPREFIX,
                                                replic_timing_file_prefix = REPCORRECTPREFIX,
                                                chrom_names = chrom_names)


get_gcCorrected_logr <- function(tumor_wgs_aliquot_id) {
  print(paste0("running sample ", tumor_wgs_aliquot_id))
  
  options(bitmapType = "cairo")
  
  TUMOURNAME <- tumor_wgs_aliquot_id

  tumourdir <- file.path(BASEOUT, paste0(TUMOURNAME, "_allelecounts"))
  
  # only continue when both alleleFrequency folders exist
  
  if (file.exists(tumourdir)) {
    setwd(tumourdir)
  } else {
    return(NULL)
  }
  
  post_correction_file <- paste0(TUMOURNAME, "_GCwindowCorrelations_afterCorrection.txt")
  if (file.exists(post_correction_file)) {
    post_correction <- read.delim(file = post_correction_file, as.is = T)
    if (any(post_correction$windowsize == "replication")) {
      return(NULL)
    }
  } else {
    return(NULL)
  }
    
  # Perform GC correction
  gc.replic.correct.wgs(Tumour_LogR_file=paste(TUMOURNAME,"_mutantLogR.tab", sep=""),
                        outfile=paste(TUMOURNAME,"_mutantLogR_gcCorrected.tab", sep=""),
                        correlations_outfile=paste(TUMOURNAME, "_GCwindowCorrelations.txt", sep=""),
                        gc_content_file_prefix=NULL,
                        replic_timing_file_prefix=NULL,
                        recalc_corr_afterwards = T,
                        chrom_names=chrom_names, 
                        gc_content_df=correctiondata$gc,
                        replic_timing_df=correctiondata$rep)
  
  return(NULL)
}



# debug(get_baf_gcCorrected_logr)

# get_baf_gcCorrected_logr(sampleid = "00c27940-c623-11e3-bf01-24c6515278c0")

gccorr <- slurm_apply(f = get_gcCorrected_logr, params = releasetable[,"tumor_wgs_aliquot_id", drop = F], jobname = "GCREcorr2", nodes = 14, cpus_per_node = 1, add_objects = ls(),
                      pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(), submit = T)
# print_job_status(gccorr)
# cancel_slurm(gccorr)

