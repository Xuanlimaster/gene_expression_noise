## Load R packages quietly
suppressPackageStartupMessages({
  library(tidyverse)
})

meta_liver <- fread(file.path("human_adult_liver", "E-MTAB-10553.cell_metadata.tsv")) %>%
  .[, 1:9] %>%
  mutate(BioSD_SAMPLE = str_extract(id, "^[^-]+")) %>%
  left_join(
    read.delim(file.path("human_adult_liver", "E-MTAB-10553.sdrf.txt")) %>%
      dplyr::select(
        BioSD_SAMPLE = "Comment.BioSD_SAMPLE.",
        BatchInfo = "Comment.ENA_SAMPLE."  # Comment.ENA_SAMPLE. as BatchInfo
      ) %>%
      distinct(BioSD_SAMPLE, BatchInfo, .keep_all = FALSE), 
    by = "BioSD_SAMPLE"
  ) %>%
  column_to_rownames("id")                 # Set cell IDs as rownames
saveRDS(meta_liver, file.path("human_adult_liver", "integrated_cell_metadata.rds"))

meta_kidney <- fread(file.path("human_adult_kidney", "E-CURD-119.cell_metadata.tsv")) %>%
  .[, 1:12] %>%
  mutate(BioSD_SAMPLE = str_extract(id, "^[^-]+")) %>%
  left_join(
    read.delim(file.path("human_adult_kidney", "E-CURD-119.sdrf.txt")) %>%
      dplyr::select(
        BioSD_SAMPLE = "Comment.BioSD_SAMPLE.",
        BatchInfo = "Comment.ENA_SAMPLE."  # Comment.ENA_SAMPLE. as BatchInfo
      ) %>%
      distinct(BioSD_SAMPLE, BatchInfo, .keep_all = FALSE), 
    by = "BioSD_SAMPLE"
  ) %>%
  column_to_rownames("id")                 # Set cell IDs as rownames
saveRDS(meta_kidney, file.path("human_adult_kidney", "integrated_cell_metadata.rds"))

meta_lung <- fread(file.path("human_adult_lung", "E-CURD-126.cell_metadata.tsv")) %>% 
  .[, 1:7] %>%
  mutate(BioSD_SAMPLE = str_extract(id, "^[^-]+")) %>%
  left_join(
    read.delim(file.path("human_adult_lung", "E-CURD-126.sdrf.txt")) %>%
      dplyr::select(
        BioSD_SAMPLE = "Comment.BioSD_SAMPLE.",
        BatchInfo = "Comment..ENA_SAMPLE." # Comment..ENA_SAMPLE. as BatchInfo
      ) %>%
      distinct(BioSD_SAMPLE, BatchInfo, .keep_all = FALSE), 
    by = "BioSD_SAMPLE"
  ) %>%
  column_to_rownames("id")                 # Set cell IDs as rownames
saveRDS(meta_lung, file.path("human_adult_lung", "integrated_cell_metadata.rds"))