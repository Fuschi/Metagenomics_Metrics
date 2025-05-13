library(tidyverse)

EU <- read_tsv("genomic_family_fragmentCountAln_ncbi_ids_EU.txt")
CPH <- read_tsv("genomic_family_fragmentCountAln_ncbi_ids_CPH.txt")
meta <- read.table("meta.tsv", header = TRUE, row.names = 1)

df_EU <- EU %>%
  separate("sample,nan", c("species","id"), sep = ",") %>%
  select(-id) %>%
  pivot_longer(-species, names_to = "sample_id", values_to = "abun") %>%
  group_by(sample_id, species) %>%
  summarise(abun = sum(abun), .groups = "drop") %>%
  pivot_wider(names_from = "sample_id", values_from = "abun") %>%
  column_to_rownames("species") %>% t()

df_CPH <- CPH %>%
  separate("sample,nan", c("species","id"), sep = ",") %>%
  select(-id) %>%
  pivot_longer(-species, names_to = "sample_id", values_to = "abun") %>%
  group_by(sample_id, species) %>%
  summarise(abun = sum(abun), .groups = "drop") %>%
  pivot_wider(names_from = "sample_id", values_from = "abun") %>%
  column_to_rownames("species") %>% t()

df <- bind_rows(as.data.frame(df_EU), as.data.frame(df_CPH))
df[is.na(df)] <- 0


idx <- purrr::map_int(rownames(meta), \(x){
  which(str_detect(rownames(df), x))
})

df_sub <- df[idx, ]
df_sub <- as_tibble(df_sub, rownames = "sample_id")

df_sub$sample_id <- str_split_fixed(df_sub$sample_id, n = 5, pattern = "_")[,3:4] %>%
  apply(1, \(x) paste(x[1], x[2], sep = "_"))


meta <- as_tibble(meta, rownames = "sample_id")
df_sub <- df_sub[rownames(meta), ]

write_tsv(df_sub, "../E-WADES_M2/counts.tsv")
write_tsv(meta, "../E-WADES_M2/meta.tsv")
