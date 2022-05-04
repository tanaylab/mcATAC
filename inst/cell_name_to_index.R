#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

cell_names <- read.table(args[1], header = FALSE, sep = "\t")[, 1]

# read from stdin
df <- read.table(file("stdin"), header = FALSE, sep = " ", col.names = c("count", "pos", "cell_name"))

# convert cell names to indices
df[, 3] <- as.numeric(factor(x = df[, 3], levels = cell_names))
df[, 2] <- as.numeric(df[, 2])

# remove non-existant cells
df <- df[!is.na(df[, 3]), ]
df <- df[, c(2, 3, 1)]
df <- dplyr::arrange(df, pos, cell_name)
# write to stdout (space separated)
write.table(
    x = df,
    file = "",
    quote = FALSE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
)
