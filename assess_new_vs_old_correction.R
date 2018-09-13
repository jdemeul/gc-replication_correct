## load all new and old pre-post data

BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/data/battenberg_logR/"

gccorrfiles <- list.files(path = BASEOUT, pattern = "*OldCorrection.txt",
                          full.names = T, recursive = T)

gccor_oldafter <- lapply(gccorrfiles, read.delim, as.is = T)
gccorrnewidxs <- which(sapply(gccor_oldafter, nrow) == 33)

gccorrfiles <- gccorrfiles[gccorrnewidxs]

gccor_oldafter <- lapply(gccorrfiles, read.delim, as.is = T)
gccor_oldafter <- do.call(rbind, lapply(gccorrfiles, read.delim, as.is = T))
gccor_oldafter$sample <- rep(sub(pattern = "_GCwindowCorrelations_afterOldCorrection.txt", replacement = "", x = basename(gccorrfiles)), each = 33)
gccor_oldafter$windowsize <- factor(gccor_oldafter$windowsize, levels = gccor_oldafter$windowsize[1:33])



# gccorr_wide <- reshape(gccor_oldafter, idvar = "sample", timevar = "windowsize", direction = "wide")


gccorrfiles_new <- sub(pattern = "Old", replacement = "", x = gccorrfiles)
gccor_newafter <- do.call(rbind, lapply(gccorrfiles_new, read.delim, as.is = T))
# gccor_newafter$sample <- rep(sub(pattern = "_GCwindowCorrelations_afterOldCorrection.txt", replacement = "", x = basename(gccorrfiles)), each = 19)
# gccor_newafter$windowsize <- factor(gccor_newafter$windowsize, levels = c(paste0(c(25,50,100,200,500), "bp"),
#                                                                           paste0(c(1,2,5,10,20,50,100,200,500), "kb"),
#                                                                           paste0(c(1,2,5,10), "Mb"),
#                                                                           "replication"))
# gccorr_wide <- cbind(gccorr_wide, reshape(gccor_newafter, idvar = "sample", timevar = "windowsize", direction = "wide"))
# # colnames(gccorr_wide)
# gccorr_wide <- gccorr_wide[, -21]
# colnames(gccorr_wide) <- c("sample", paste0("old.", colnames(gccorr_wide)[2:20]), paste0("new.", colnames(gccorr_wide)[21:39]))
# colnames(gccorr_wide) <- sub(pattern = ".1$", replacement = "", x = colnames(gccorr_wide))

# write.table(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/GCwindowCorrelations_acrossPCAWG_withrep.txt",
#             x = oldgccorr, sep = "\t", quote = F, row.names = F)

gccorr_diff <- gccor_oldafter
gccorr_diff$correlation <- gccor_oldafter$correlation - gccor_newafter$correlation

library(ggplot2)
p1 <- ggplot(data = gccorr_diff, mapping = aes(x = windowsize, y = correlation, group = sample)) 
p1 <- p1 + geom_line(alpha = .05) + theme_minimal()
p1

p1 <- ggplot(data = gccorr_diff, mapping = aes(x = windowsize, y = correlation)) 
p1 <- p1 + geom_boxplot() + theme_minimal() + labs(y = "Correlation post old - post new correction", title = "Results on 250 PCAWG samples")
p1
ggsave(filename = "boxplot_corr_post_old-post_new_correct.png", plot = p1)
