# Step 15: the mammal tree of the student project: the dated tree of Upham et al. (2019; DNA-only, node-dated,
# maximum clade credibility) cut down to the species of the student's files, with the tips named as in them.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "lib.R"))
suppressMessages({library(ape); library(jsonlite)})
tr <- load_tree()
m <- fromJSON(file.path(OUT, "work", "tree_tip_map.json"))  # tree binomial (Genus_species) -> name in the data
names(m) <- gsub("_", " ", names(m))
tips <- intersect(names(m), tr$tip.label)
phy <- keep.tip(tr, tips)
phy$tip.label <- unname(unlist(m[phy$tip.label]))
phy$node.label <- NULL
phy <- ladderize(phy)
dir.create(file.path(OUT, "student", "fnss", "tree"), showWarnings = FALSE, recursive = TRUE)
write.tree(phy, file.path(OUT, "student", "fnss", "tree", "mammals.nwk"), digits = 4)
cat(length(phy$tip.label), "species in the tree;", length(setdiff(names(m), tr$tip.label)), "not in it:",
    paste(head(setdiff(names(m), tr$tip.label), 30), collapse = ", "), "\n")
