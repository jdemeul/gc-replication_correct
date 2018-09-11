## check possibility of large-scale windows being affected by CNAs
library(ggplot2)

# go for ~ approx diploid pure tumour and lose/gain GC-poor chrom
purity_ploidy_file <- "/srv/shared/vanloo/ICGC_consensus_copynumber/consensus.20170217.purity.ploidy.txt.gz"
purity_ploidy <- read.delim(file = purity_ploidy_file, as.is = T)

# f8f7f274-dd98-3bb0-e040-11ac0c483fcd is highly pure, diploid sample
# check effect of a loss on correlations

BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/data/battenberg_logR/"
chrom_names <- c(1:22, "X")

TUMOURNAME <- "f8f7f274-dd98-3bb0-e040-11ac0c483fcd"
tumourdir <- file.path(BASEOUT, paste0(TUMOURNAME, "_allelecounts"))
setwd(tumourdir)

gc_sub <- correctiondata$gc[seq(from = 1, to = nrow(correctiondata$gc), length.out = 100000),
                            c("chr", "Position", "10Mb")]

plot(x = 1:10000, y = correctiondata$gc$`10Mb`[])

p1 <- ggplot(data = gc_sub, mapping = aes(x = Position, y = `10Mb`)) + geom_point() + facet_wrap(~chr)
p1

# GC poor chroms: (q-arm of) 4, 5, 13
# GC rich chroms: (p-arm of) 1, 19, 17, 16

sort(c(by(data = gc_sub$`10Mb`, INDICES = gc_sub$chr, FUN = mean, na.rm =T)))


### do logr correction
Tumour_LogR_file=paste(TUMOURNAME,"_mutantLogR.tab", sep="")
outfile=paste(TUMOURNAME,"_mutantLogR_gcCorrected.tab", sep="")
correlations_outfile=paste(TUMOURNAME, "_GCwindowCorrelations.txt", sep="")
gc_content_file_prefix=NULL
replic_timing_file_prefix=NULL
recalc_corr_afterwards = T
chrom_names=chrom_names
gc_content_df=correctiondata$gc
replic_timing_df=correctiondata$rep



Tumor_LogR = readr::read_tsv(file = Tumour_LogR_file, col_names = T, col_types = "cin")

GC_data <- gc_content_df
rm(gc_content_df)

replic_data <- replic_timing_df
rm(replic_timing_df)

# omit non-matching loci, replication data generated at exactly same GC loci
locimatches <- match(x = paste0(Tumor_LogR$Chromosome, "_", Tumor_LogR$Position),
                     table = paste0(GC_data$chr, "_", GC_data$Position))
Tumor_LogR <- Tumor_LogR[which(!is.na(locimatches)), ]
GC_data <- GC_data[na.omit(locimatches), ]
replic_data <- replic_data[na.omit(locimatches), ]
rm(locimatches)

corr = abs(cor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[,3], use="complete.obs")[,1])
corr_rep = abs(cor(replic_data$timing, Tumor_LogR[,3], use="complete.obs")[,1])

index_1kb = which(names(corr)=="1kb")
maxGCcol_insert = names(which.max(corr[1:index_1kb]))

index_100kb = which(names(corr)=="100kb")
# start large window sizes at 5kb rather than 2kb to avoid overly correlated expl variables
maxGCcol_amplic = names(which.max(corr[(index_1kb+2):index_100kb]))

# Multiple regression 
corrdata <- data.frame(logr = Tumor_LogR[,3, drop = T],
                       GC_insert = GC_data[,maxGCcol_insert, drop = T],
                       GC_amplic = GC_data[,maxGCcol_amplic, drop = T],
                       replic = replic_data[, "timing", drop = T])
model = lm(logr ~ poly(GC_insert, 2, raw = T) + poly(GC_amplic, 2, raw = T) + poly(replic, 2, raw = T), data = corrdata, na.action="na.exclude")
model = lm(logr ~ ns(x = GC_insert, df = 5, intercept = T) + ns(x = GC_amplic, df = 5, intercept = T) + ns(x = GC_amplic, df = 5, intercept = T), data = corrdata, na.action="na.exclude")

Tumor_LogR_corrected = residuals(model)
rm(model, corrdata)

corr = data.frame(windowsize=c(names(corr), "replication"), correlation=c(corr, corr_rep))
corr


# Recalculate the correlations to see how much there is left
corr = abs(cor(GC_data[, 3:ncol(GC_data)], Tumor_LogR_corrected, use="complete.obs")[,1])
corr_rep = abs(cor(replic_data$timing, Tumor_LogR_corrected, use="complete.obs"))
corr = data.frame(windowsize=c(names(corr), "replication"), correlation=c(corr, corr_rep))
corr
# write.table(corr, file=gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)



###############################
## fake loss of chrom 4 q-arm
Tumor_LogR_mod <- Tumor_LogR
Tumor_LogR_mod$`f8f7f274-dd98-3bb0-e040-11ac0c483fcd` <- ifelse(Tumor_LogR_mod$Chromosome == "4" & Tumor_LogR_mod$Position > 5e7, Tumor_LogR_mod$`f8f7f274-dd98-3bb0-e040-11ac0c483fcd` - .97, Tumor_LogR_mod$`f8f7f274-dd98-3bb0-e040-11ac0c483fcd`)

corr = abs(cor(GC_data[, 3:ncol(GC_data)], Tumor_LogR_mod[,3], use="complete.obs")[,1])
corr_rep = abs(cor(replic_data$timing, Tumor_LogR_mod[,3], use="complete.obs")[,1])

corr_orig = abs(cor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[,3], use="complete.obs")[,1])
corr_rep_orig = abs(cor(replic_data$timing, Tumor_LogR[,3], use="complete.obs")[,1])

corr = data.frame(windowsize=c(names(corr), "replication"), correlation=c(corr, corr_rep))
corr_orig = data.frame(windowsize=c(names(corr_orig), "replication"), correlation=c(corr_orig, corr_rep_orig))

corr$truecorr <- corr_orig$correlation
corr$windowsize <- factor(corr$windowsize, levels = corr$windowsize)

p1 <- ggplot(data = corr) + geom_point(mapping = aes(x = windowsize, y = correlation), colour = "red") + geom_point(mapping = aes(x = windowsize, y = truecorr), colour = "blue") + ylim(c(0,.25))
p1
