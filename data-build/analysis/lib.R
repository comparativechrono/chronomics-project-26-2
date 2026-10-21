# Shared R helpers for the FNSS analysis
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
ROOT <- if (length(script_path)) normalizePath(file.path(dirname(script_path), "..")) else normalizePath(".")
DATA <- file.path(ROOT, "data")
OUT <- file.path(ROOT, "outputs")
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

rtsv <- function(path, ...) read.delim(path, sep = "\t", quote = "", comment.char = "",
                                        stringsAsFactors = FALSE, na.strings = c("", "NA"), check.names = FALSE, ...)
wtsv <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# the dated mammal tree (Upham et al. 2019, DNA-only, node-dated, MCC) with tips named "Genus species"
load_tree <- function() {
  suppressMessages(library(ape))
  tr <- read.nexus(file.path(DATA, "tree", "MamPhy_DNAonly_4098sp_NDexp_MCC.tre"))
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  tr <- drop.tip(tr, grep("^_", tr$tip.label, value = TRUE))  # outgroup
  p <- strsplit(tr$tip.label, "_")
  tr$tip.label <- vapply(p, function(x) paste(x[1], x[2]), "")
  # the tip labels of the source tree carry family and order (Genus_species_FAMILY_ORDER)
  attr(tr, "order_of") <- setNames(vapply(p, function(x) x[length(x)], ""), tr$tip.label)
  tr
}

# names in the data -> tree binomial
tree_name <- function(sp) {
  syn <- rtsv(file.path(DATA, "tree", "tree_name_synonyms.tsv"))
  m <- match(sp, syn$name)
  out <- sp
  out[!is.na(m)] <- gsub("_", " ", syn$tree_binomial[m[!is.na(m)]])
  out
}

fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3)))
