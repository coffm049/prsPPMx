options(tidyverse.quiet = TRUE, readr.show_col_types = FALSE)
library(tidyverse)

files <- list.files("simulation/results", pattern = "confound", full.names = TRUE)

read_results <- function(filepath, csv = T) {
  { if (csv) {
    read_csv(filepath, col_names = c("f1", "R2", "confounding"), show_col_types = F)
  } else {
    read_table(filepath, col_names = c("f1", "R2", "confounding"), show_col_types = F)
  }} %>%
    mutate(f1 = as.numeric(f1), R2 = as.numeric(R2))
}

df <- files %>%
  map(function(x) {
    read_results(x, csv = T) %>%
      mutate(., name = x)
  }, .progress = TRUE) %>%
  reduce(bind_rows) %>%
  drop_na() %>%
  mutate(
    name = stringr::str_sub(basename(name), start = 9, end = 11),
    name = ifelse(name == "0.2", "0.25", name)
  )

df %>%
  ggplot(aes(x = name, y = f1, fill = confounding)) +
  geom_boxplot() +
  ggtitle("F1 Score Comparison") +
  xlab("Projection") +
  ylab("F1 Score")

ggsave("f1_comparison.png", width = 5, height = 5, dpi = 300)
