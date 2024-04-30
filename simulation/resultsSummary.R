options(tidyverse.quiet = TRUE, readr.show_col_types = FALSE)
library(tidyverse)
library(janitor)

df <- read_table("temp/simulation.results") %>%
    clean_names() %>%
    filter(!grepl("f1", f1)) %>%
    mutate(across(c(f1, r2), as.numeric)) 


df %>%
    ggplot(aes(x = proj, y = f1)) +
    geom_boxplot()

ggsave("temp/f1_comparison.png", width = 5, height = 5, dpi = 300)
