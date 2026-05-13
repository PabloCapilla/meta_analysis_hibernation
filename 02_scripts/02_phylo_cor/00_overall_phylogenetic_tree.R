###
###
#' Script to prepare an overall phylogenetic correlation matrix
#' using all species present in the processed meta-analysis datasets
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
library(ggtree)
library(ggimage)
library(treeio)


##
##### data #####
##
data_folder <- "./01_data/processed_RDS_data_files"

data_files <- list.files(
  path = data_folder,
  pattern = "_metaanalysis_.*_full_data\\.RDS$",
  full.names = TRUE
)

if (length(data_files) == 0) {
  stop(
    "No meta-analysis full-data RDS files were found in ./01_data/processed_RDS_data_files"
  )
}

data_list <- lapply(data_files, readRDS)

missing_species_column <- data_files[
  !vapply(
    data_list,
    function(x) is.data.frame(x) && "Species" %in% names(x),
    logical(1)
  )
]

if (length(missing_species_column) > 0) {
  stop(
    "These files do not contain a Species column: ",
    paste(basename(missing_species_column), collapse = ", ")
  )
}

missing_ref_column <- data_files[
  !vapply(
    data_list,
    function(x) is.data.frame(x) && "Reference" %in% names(x),
    logical(1)
  )
]

if (length(missing_ref_column) > 0) {
  stop(
    "These files do not contain a Reference column: ",
    paste(basename(missing_ref_column), collapse = ", ")
  )
}

data <- bind_rows(data_list, .id = "source_file") %>%
  mutate(
    source_file = basename(data_files[as.integer(source_file)])
  ) |>
  select(StudyID, Species, Reference, source_file)

## number of studies
length(unique(data$Reference))
length(unique(data$Species))

data |>
  filter(is.na(Reference))

nrow(data)
head(data)
str(data)

list_studies <- data |>
  select(Reference) |>
  group_by(Reference) |>
  filter(row_number() == 1) |>
  arrange(Reference)

write.csv(list_studies, './01_data/list_unique_studies.csv')

##
##### recovering phylogenetic relationship #####
##
data$Species <- replace_values(
  data$Species,
  "Myotis_pilosus" ~ "Myotis_ricketti",
  "Lithobates pipiens" ~ "Rana pipiens",
  "Ranoidea alboguttata" ~ "Cyclorana alboguttata" # correcting synomyn names manuallty
)

species_data <- data %>%
  filter(!is.na(Species), Species != "") %>%
  distinct(Species) %>%
  mutate(
    scientific_name = gsub(
      x = Species,
      pattern = "_",
      replacement = " "
    )
  )

taxa.corrected <- tnrs_match_names(names = species_data$scientific_name)


# check approximate matches
taxa.corrected[taxa.corrected$approximate_match == TRUE, ]

# check synonym matches
taxa.corrected[taxa.corrected$is_synonym == TRUE, ]

# check number of matches
taxa.corrected[taxa.corrected$number_matches > 1, ]

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
tree <- tol_induced_subtree(
  ott_ids = taxa.corrected[["ott_id"]],
  label_format = "name"
)

tree$tip.label[tree$tip.label == "mrcaott307132ott864246"] <- "Myotis lucifugus"


##
##
##### Tree plot and proportions of data per species #####
##
##

overall_tree <- ggtree(tree, size = 1.2) +
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

overall_tree

ggsave(
  filename = here(
    "./03_plots/overall_phylo_tree.png"
  ),
  plot = overall_tree,
  device = "png",
  height = 200,
  width = 125,
  units = "mm"
)

ggsave(
  filename = here(
    "./03_plots/overall_phylo_tree.pdf"
  ),
  plot = overall_tree,
  height = 200,
  width = 125,
  units = "mm"
)

#####

##
##### list of all papers included #####
##
data |>
  distinct(Reference)
