# check GC correction correlation windows
library(ggplot2)
# library(ggforce)

BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/data/battenberg_logR/"

gccorrfiles <- list.files(path = BASEOUT, pattern = "*_GCwindowCorrelations_beforeCorrection.txt",
                          full.names = T, recursive = T)

gccor <- do.call(rbind, lapply(gccorrfiles, read.delim, as.is = T))
gccor$sample <- rep(sub(pattern = "_GCwindowCorrelations_beforeCorrection.txt", replacement = "", x = basename(gccorrfiles)), each = 18)
gccor$windowsize <- factor(gccor$windowsize, levels = c(paste0("X", c(25,50,100,200,500), "bp"),
                                                        paste0("X", c(1,2,5,10,20,50,100,200,500), "kb"),
                                                        paste0("X", c(1,2,5,10), "Mb")))

gccorr_wide <- reshape(gccor, idvar = "sample", timevar = "windowsize", direction = "wide")
write.table(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/GCwindowCorrelations_acrossPCAWG.txt",
            x = gccorr_wide, sep = "\t", quote = F, row.names = F)

gccorr_wide$max_short <- levels(gccor$windowsize)[1:14][c(apply(X = gccorr_wide[, 2:15], MARGIN = 1, which.max))]
gccorr_wide$max_long <- levels(gccor$windowsize)[15:18][c(apply(X = gccorr_wide[, 16:19], MARGIN = 1, which.max))]

gccorr_wide$max_short <- factor(gccorr_wide$max_short, levels = c(paste0("X", c(25,50,100,200,500), "bp"),
                                                            paste0("X", c(1,2,5,10,20,50,100,200,500), "kb")))
gccorr_wide$max_long <- factor(gccorr_wide$max_long, levels = paste0("X", c(1,2,5,10), "Mb"))



table(gccorr_wide$max_short)
table(gccorr_wide$max_long)



p1 <- ggplot(data = gccorr_wide[gccorr_wide$max_short == "X500bp", ], mapping = aes(x = correlation.X500bp, y = correlation.X1Mb)) 
p1 <- p1 + geom_point(alpha = .05) + theme_minimal()
p1

p1 <- ggplot(data = gccorr_wide, mapping = aes(x = gccorr_wide$max_short)) 
p1 <- p1 + geom_bar() + theme_minimal()
p1

p1 <- ggplot(data = gccorr_wide, mapping = aes(x = gccorr_wide$max_long)) 
p1 <- p1 + geom_bar() + theme_minimal()
p1


p1 <- ggplot(data = gccor, mapping = aes(x = windowsize, y = correlation, group = sample)) 
p1 <- p1 + geom_line(alpha = .05) + theme_minimal()
p1

p1 <- ggplot(data = gccor, mapping = aes(x = windowsize, y = correlation)) 
p1 <- p1 + geom_jitter(alpha = .1, width = .4, height = 0) + scale_y_log10() + geom_boxplot(alpha = .5, outlier.shape = NA)
p1

p1 <- ggplot(data = gccor, mapping = aes(x = windowsize, y = correlation)) + geom_violin(alpha = 1, draw_quantiles = c(.5))
p1 <- p1 + geom_jitter(alpha = .05, width = .35, height = 0, size = 1) + scale_y_log10()
p1 <- p1 + theme_minimal()
p1

temp <- c(by(data = gccor, INDICES = gccor$sample, FUN = function(x) levels(x$windowsize)[which.max(x$correlation)]))
hist(temp)

hist(c(by(data = gccor, INDICES = gccor$sample, FUN = function(x) x[which.max(x$correlation), "windowsize"])), breaks = 0:18, labels = levels(gccor$windowsize), xlab = "", main = "Histogram of max correlation")

hist(c(by(data = gccor, INDICES = gccor$sample, FUN = function(x) mean(x$correlation))), breaks = 50,xlab = "", main = "Histogram of max correlation")


