library(tidyverse)

process_strata <- function(df) {
  read_tsv(df, col_types = cols(.default = "c")) %>% 
    separate(INFO, into = c("AC", "AN"), sep = ";") %>% 
    select(-AN, -ID, -QUAL, -FILTER) %>% 
    mutate(across(AC, str_sub, start = 4L)) %>%
    separate_rows(ALT, AC, sep = ",") %>% 
    distinct(across(`#CHROM`:ALT), .keep_all = TRUE)
}

snakemake@input %>% 
  map(process_strata) %>% 
  set_names(c("ALL", "FEMALE", "MALE")) %>% 
  bind_rows(.id = "STRATA") %>% 
  pivot_wider(names_from = STRATA, values_from = AC) %>% 
  write_tsv(snakemake@output[[1]]) 
