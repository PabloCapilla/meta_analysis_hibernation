###
###
#' Script to prepare phylogenetic correlation matrix
#'
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list = ls())


##
##### libraries #####
##
library(dplyr)
library(tidyr)
library(metafor)
library(rotl)
library(ape)
library(phytools)

#if (!require("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}

#BiocManager::install("ggtree")
library(ggtree)

#BiocManager::install("ggimage", force = T)
library(ggimage)

#BiocManager::install("treeio", force = T)
library(treeio)

##
##### data #####
##
data <- readRDS(
  "./01_data/processed_RDS_data_files/06_GC_metaanalysis_arousal_full_data.RDS"
)
head(data)
str(data)

##
##### recovering phylogenetic relationship #####
##
data$scientific_name <- gsub(
  x = data$Species,
  pattern = "_",
  replacement = " "
)

taxa.corrected <- tnrs_match_names(names = unique(data$scientific_name))


# check approximate matches: OK
taxa.corrected[taxa.corrected$approximate_match == TRUE, ]

# check synonyms matches: OK
taxa.corrected[taxa.corrected$is_synonym == TRUE, ]

# check number of matches: OK
taxa.corrected[taxa.corrected$number_matches > 1, ]

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
tree <- tol_induced_subtree(
  ott_ids = taxa.corrected[["ott_id"]],
  label_format = "name"
)

tree$tip.label <- replace_values(
  x = tree$tip.label,
  "mrcaott307132ott864246" ~ "Myotis lucifugus"
)

# save tree
saveRDS(
  tree,
  file = "./01_data/processed_RDS_data_files/06_GC_phylogenetic_tree_arousal.RDS"
)
##
##
##### Tree plot and proportions of data per species #####
##
##

p2 <- ggtree(tree, size = 1.2) +
  xlim_tree(c(0, 35)) +
  geom_tiplab(
    aes(
      label = paste0(
        "italic('",
        sub(x = label, pattern = "_", replacement = " "),
        "')"
      )
    ),
    parse = TRUE,
    offset = 0.5,
    size = 5,
    hjust = 0
  )

p2
##
##
##### Creating phylogenetic relationship matrix for meta-analyses #####
##
##

# check for the existence of polytomies
is.binary.tree(tree) # there are no polytomies


# checking that the tree includes every species in data table
tree$tip.label <- gsub("_", " ", tree$tip.label)
intersect(as.character(tree$tip.label), as.character(data$scientific_name))
setdiff(as.character(data$scientific_name), as.character(tree$tip.label)) #listed in our database but not in the tree
setdiff(as.character(tree$tip.label), as.character(data$scientific_name)) # listed in the tree but not in our database

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# matrix to be included in the models
phylo_cor <- vcv(phylo_branch, cor = T)

#
# save matrix for analyses
saveRDS(
  phylo_cor,
  file = "./01_data/processed_RDS_data_files/06_GC_phylogenetic_correlations_arousal.RDS"
)
