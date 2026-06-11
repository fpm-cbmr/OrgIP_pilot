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


# Heatmap of markers ------------------------------------------------------

# Myoblasts:

heatmap_mb <- data_long |>
  filter(
    diff == "mb"
  ) |>
  filter(protein %in% c("Lamp1",
                        "Tomm20",
                        "Atp2a2",
                        "Tfrc")) |>
  group_by(protein, tag) |>
  summarise(mean_LFQ = mean(LFQ, na.rm = T)) |>
  ungroup() |>
  tidyr::pivot_wider(names_from = protein,
                     values_from = mean_LFQ) |>
  tibble::column_to_rownames("tag") |>
  scale() |>
  as.data.frame() |>
  tibble::rownames_to_column("tag") |>
  tidyr::pivot_longer(
    cols = !tag,
    names_to = "protein",
    values_to = "scaled_LFQs"
  )

heatmap_mb |>
  ggplot(
    aes(
      x = tag,
      y = protein,
      fill = scaled_LFQs,
    )
  ) +
  ggtitle("LFQ intensity of IP baits - Myoblasts") +
  geom_tile(color = "black", alpha = 0.85) +
  scale_fill_viridis_c(option = "plasma") +
  theme_classic() +
  ylab("MS-quantified protein") +
  xlab("Protein tag for IP") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    text = element_text(size = 8),
    legend.key.size = unit(2, "mm")
  )

# ggsave(here::here("plots/heatmap_ip_baits_mb.pdf"),
#        units = "mm",
#        height = 60,
#        width = 120)

# The IP for Tomm20 and Lamp1 in myoblasts works

# Myotubes:

heatmap_mt <- data_long |>
  filter(
    diff == "mt"
  ) |>
  filter(protein %in% c("Lamp1",
                        "Tomm20",
                        "Atp2a2",
                        "Tfrc")) |>
  group_by(protein, tag) |>
  summarise(mean_LFQ = mean(LFQ, na.rm = T)) |>
  ungroup() |>
  tidyr::pivot_wider(names_from = protein,
                     values_from = mean_LFQ) |>
  tibble::column_to_rownames("tag") |>
  scale() |>
  as.data.frame() |>
  tibble::rownames_to_column("tag") |>
  tidyr::pivot_longer(
    cols = !tag,
    names_to = "protein",
    values_to = "scaled_LFQs"
  )

heatmap_mt |>
  ggplot(
    aes(
      x = tag,
      y = protein,
      fill = scaled_LFQs,
    )
  ) +
  ggtitle("LFQ intensity of IP baits - Myotubes") +
  geom_tile(color = "black", alpha = 0.85) +
  scale_fill_viridis_c(option = "plasma") +
  theme_classic() +
  ylab("MS-quantified protein") +
  xlab("Protein tag for IP") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    text = element_text(size = 8),
    legend.key.size = unit(2, "mm")
  )

# ggsave(here::here("plots/heatmap_ip_baits_mt.pdf"),
#        units = "mm",
#        height = 60,
#        width = 120)


# Organelle enrichment with IPs -------------------------------------------

# Load cell localization annotations:

annotations <- vroom::vroom(here::here("data/Hein_cell_localization_annotations.csv")) |>
  mutate(Gene_name_canonical = stringr::str_to_title(Gene_name_canonical)) |>
  mutate(protein = Gene_name_canonical)

data_annotated <- data_long |>
  inner_join(annotations)

# Mitochondria:

# Myoblasts:

data_mito <- data_annotated |>
  dplyr::filter(
graph_localization_annotation == "mitochondrion"
  ) |>
  filter(diff == "mb") |>
  group_by(protein, tag) |>
  summarise(mean_LFQ = mean(LFQ, na.rm = T)) |>
  ungroup()

data_mito |>
  ggplot(
    aes(
      x = tag,
      y = mean_LFQ,
      fill = tag
    )
  ) +
  geom_violin(alpha = 0.75) +
  geom_point(size = 0.5, position = position_jitter(width = 0.1), size = 1) +
  geom_point(
    data = ~filter(.x, protein == "Tomm20"),
    color = "red", size = 3, shape = 21, stroke = 1.5,
    show.legend = FALSE
  ) +
  ggtitle("LFQ intensity of Mitochondrial proteins - Myoblasts") +
  theme_classic() +
  scale_fill_viridis_d(option = "turbo") +
  labs(caption = "Red dot: Tomm20") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    text = element_text(size = 8),
    legend.key.size = unit(3, "mm")
  )

ggsave(here::here("plots/violins_mito_mb.pdf"),
       units = "mm",
       height = 60,
       width = 120)


# Myotubes:

data_mito <- data_annotated |>
  dplyr::filter(
    graph_localization_annotation == "mitochondrion"
  ) |>
  filter(diff == "mt") |>
  group_by(protein, tag) |>
  summarise(mean_LFQ = mean(LFQ, na.rm = T)) |>
  ungroup()

data_mito |>
  ggplot(
    aes(
      x = tag,
      y = mean_LFQ,
      fill = tag
    )
  ) +
  geom_violin(alpha = 0.75) +
  geom_point(size = 0.5, position = position_jitter(width = 0.1), size = 1) +
  geom_point(
    data = ~filter(.x, protein == "Tomm20"),
    color = "red", size = 3, shape = 21, stroke = 1.5,
    show.legend = FALSE
  ) +
  ggtitle("LFQ intensity of Mitochondrial proteins - Myotubes") +
  theme_classic() +
  scale_fill_viridis_d(option = "turbo") +
  labs(caption = "Red dot: Tomm20") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    text = element_text(size = 8),
    legend.key.size = unit(3, "mm")
  )

ggsave(here::here("plots/violins_mito_mt.pdf"),
       units = "mm",
       height = 60,
       width = 120)

# Lysosomes:

# Myoblasts:

data_lyso <- data_annotated |>
  dplyr::filter(
    graph_localization_annotation == "lysosome"
  ) |>
  filter(diff == "mb") |>
  group_by(protein, tag) |>
  summarise(mean_LFQ = mean(LFQ, na.rm = T)) |>
  ungroup()

data_lyso |>
  ggplot(
    aes(
      x = tag,
      y = mean_LFQ,
      fill = tag
    )
  ) +
  geom_violin(alpha = 0.75) +
  geom_point(size = 0.5, position = position_jitter(width = 0.1), size = 1) +
  geom_point(
    data = ~filter(.x, protein == "Lamp1"),
    color = "red", size = 3, shape = 21, stroke = 1.5,
    show.legend = FALSE
  ) +
  ggtitle("LFQ intensity of Lysosomal proteins - Myoblasts") +
  theme_classic() +
  scale_fill_viridis_d(option = "turbo") +
  labs(caption = "Red dot: Tomm20") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    text = element_text(size = 8),
    legend.key.size = unit(3, "mm")
  )

ggsave(here::here("plots/violins_lyso_mb.pdf"),
       units = "mm",
       height = 60,
       width = 120)

# Myotubes:

data_lyso <- data_annotated |>
  dplyr::filter(
    graph_localization_annotation == "lysosome"
  ) |>
  filter(diff == "mt") |>
  group_by(protein, tag) |>
  summarise(mean_LFQ = mean(LFQ, na.rm = T)) |>
  ungroup()

data_lyso |>
  ggplot(
    aes(
      x = tag,
      y = mean_LFQ,
      fill = tag
    )
  ) +
  geom_violin(alpha = 0.75) +
  geom_point(size = 0.5, position = position_jitter(width = 0.1), size = 1) +
  geom_point(
    data = ~filter(.x, protein == "Lamp1"),
    color = "red", size = 3, shape = 21, stroke = 1.5,
    show.legend = FALSE
  ) +
  ggtitle("LFQ intensity of Lysosomal proteins - Myoblasts") +
  theme_classic() +
  scale_fill_viridis_d(option = "turbo") +
  labs(caption = "Red dot: Tomm20") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    text = element_text(size = 8),
    legend.key.size = unit(3, "mm")
  )

ggsave(here::here("plots/violins_lyso_mt.pdf"),
       units = "mm",
       height = 60,
       width = 120)
