library(tidyverse)

data <- vroom::vroom(here::here("data/data_named.csv")) |>
  select(!1) |>
  tibble::column_to_rownames("sample_id") |>
  log2()

data_long <- data |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(
    cols = !sample_id,
    names_to = "protein",
    values_to = "LFQ"
  ) |>
  tidyr::extract(
    col = sample_id,
    into = c("tag", "diff", "replicate"),
    regex = "^(.+)_(mb|mt)_([0-9]+)$",
    remove = FALSE
  ) |>
  mutate(
    replicate = paste0("R", replicate),
    tag = gsub("(?<=[a-z])_(?=[0-9])|(?<=[0-9])_(?=[a-z])", "", tag, perl = TRUE)
  ) |>
  dplyr::select(!sample_id)


data_long |>
  filter(protein == "Tomm20") |>
  ggplot(
    aes(
      x = tag,
      y = LFQ,
      fill = tag
    )
  ) +
  geom_boxplot() +
  facet_grid(~diff)
