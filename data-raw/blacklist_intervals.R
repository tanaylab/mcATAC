## Generate ENCODE blacklist intervals for hg38, hg19, and mm10

library(glue)
library(tidyverse)
library(misha.ext)
library(misha)
library(tgutil)

wd <- setwd(tempdir())
on.exit(setwd(wd))
system("git clone git@github.com:Boyle-Lab/Blacklist.git")

setwd("Blacklist")


create_blacklist_intervals <- function(genome) {
    lst <- fread(glue("lists/{genome}-blacklist.v2.bed.gz"), col.names = c("chrom", "start", "end", "type")) %>% as_tibble()
    gset_genome(genome)
    gdir.create("ENCODE")
    gintervals.save("ENCODE.blacklist", lst)
    gdb.reload()
}

create_blacklist_intervals("hg38")
create_blacklist_intervals("hg19")
create_blacklist_intervals("mm10")
