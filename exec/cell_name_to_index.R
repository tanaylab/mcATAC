#!/usr/bin/env Rscript

library(dplyr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

cell_names <- read.table(args[1], header = FALSE, sep = "\t")[, 1]
start <- as.numeric(args[2])
end <- as.numeric(args[3])

# read from stdin
df <- read.table(file("stdin"), header = FALSE, sep = " ", col.names = c("count", "pos", "cell_name"))

df <- df %>%
    # note that we include the last base so this works only if the start and end parameters are 0-based
    filter(pos >= start, pos <= end) %>% 
    mutate(
        # convert cell names to indices
        cell_name = as.numeric(factor(x = cell_name, levels = cell_names)),
        # adjust position to start
        pos = as.numeric(pos) - start
    ) %>%
    # remove non-existant cells and positions which are negative (they would be present in other regions)
    filter(!is.na(cell_name), pos > 0) %>%
    select(pos, cell_name, count) %>%
    arrange(pos, cell_name)

# write to stdout (space separated)
write.table(
    x = df,
    file = "",
    quote = FALSE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
)
