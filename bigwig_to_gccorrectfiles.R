### check replication timing contribution
library(readr)
library(rtracklayer)

## read wavelet smoothed repli-seq timing data
waveletfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2018_GC_correction/data/UW_repliseq_wavelet/", pattern = "*.bigWig", full.names = T)
wavelets <- lapply(waveletfiles, import.bw)

# check size of each bw
# lapply(wavelets, length)

## start on final df and average scores
avg_wave <- wavelets[[1]]
mcols(avg_wave)$score <- colMeans(do.call(what = rbind, args = lapply(X = wavelets, FUN = function(x) mcols(x)$score)))


## read loci files used in Battenberg and turn into GRanges
locifiles <- paste0("/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_loci/1000genomesloci2012_chr", 1:23, ".txt")
loci <- do.call(rbind, lapply(locifiles, FUN = read_tsv, col_types = "ci", col_names = c("chr", "pos")))
loci_gr <- GRanges(seqnames = loci$chr, ranges = IRanges(start = loci$pos, width = 1))
rm(loci)

## change chr naming and annotate SNP loci with repli-seq score 
seqlevelsStyle(avg_wave) <- "Ensembl"
mcols(loci_gr)$score <- mcols(avg_wave)[nearest(x = loci_gr, subject = avg_wave, select = "arbitrary") , "score"]

## Quick QC, chekc those 2 SNPs with negative timing score
# plot(start(loci_gr)[20950900:20951000], mcols(loci_gr)[20950900:20951000, "score"])
# plot(start(loci_gr)[20950900:20951000], mcols(loci_gr)[20950900:20951000, "score"])
# summary(mcols(loci_gr)$score[-which(mcols(loci_gr)$score < 0)])
# loci_gr[which(mcols(loci_gr)$score < 0)]
# loci_gr[20950915:20950925]


## write out correction files
loci_gr_split <- split(x = loci_gr, f = seqnames(loci_gr))
lapply(loci_gr_split, FUN = function(x) {
  write_tsv(path = paste0("/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_2012_v3_repliTiming/1000_genomes_replication_timing_chr_", seqnames(x[1]), ".txt"),
            x = data.frame(chr = seqnames(x), pos = start(x), timing = mcols(x)$score))
  return(NULL)
  })

