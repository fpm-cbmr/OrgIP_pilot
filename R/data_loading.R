# Load libraried

library(tidyverse)

raw_data <- vroom::vroom(here::here("data-raw/20260610_174505_OrgIP_pilot_myoblasts_v_myotubes_Report.tsv"))

metadata <- vroom::vroom(here::here("data/metadata.csv")) |>
  select(c(
    `Analytical Sample ID`,
    `Biological Sample ID`,
    `Plate Position`
  ))

colnames(metadata) <- snakecase::to_snake_case(colnames(metadata))


data <- raw_data |>
  select(!PG.ProteinGroups) |>
  filter(!duplicated(PG.Genes)) |>
  filter(!PG.Genes == "") |>
  tibble::column_to_rownames("PG.Genes") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("plate_position") |>
  mutate(plate_position = gsub(".*_([A-H](0?[1-9]|1[0-2]))_.*", "\\1", plate_position)) |>
  inner_join(metadata) |>
  mutate(sample_id = biological_sample_id) |>
  dplyr::select(!c(plate_position, analytical_sample_id, biological_sample_id)) |>
  dplyr::select(sample_id, everything()) |>
  dplyr::mutate(
    sample_id = snakecase::to_snake_case(sample_id)
  )

write.csv(data, here::here("data/data_named.csv"))
