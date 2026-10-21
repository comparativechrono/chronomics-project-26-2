# 08_msleep.R FILE – write msleep (Savage & West 2007; the data set of the R package ggplot2) as a tab-separated
# table, empty cells for NA.
out <- commandArgs(trailingOnly = TRUE)[1]
data("msleep", package = "ggplot2")
write.table(as.data.frame(msleep), out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(nrow(msleep), "species in msleep\n")
