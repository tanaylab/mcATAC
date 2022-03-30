
<!-- README.md is generated from README.Rmd. Please edit that file -->
# mcATAC

## Installation

You can install the development version of mcATAC like so:

``` r
remotes::install_github("tanaylab/mcATAC")
```

## Example

``` r
if (!dir.exists("pbmc_data")){
  download_pbmc_example_data()
}
```

### Import ATAC dataset

``` r
atac_sc <- import_from_10x("pbmc_data", genome = "hg38", id = "PBMC", description = "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)")
#>   • Importing matrix
#> ℹ Imported a matrix of 11909 cells and 144978 features
#>   • Importing features
#> ℹ Removed 107861 ATAC peaks which were all zero
#> ℹ 107861 ATAC peaks
#> ! removed 32 peaks from the following chromosome(s) which are missing from hg38: 'KI270727.1, GL000194.1, GL000205.2, GL000195.1, GL000219.1, KI270734.1, KI270721.1, KI270726.1, KI270713.1'
#> ✔ successfully imported to an ScATAC object with 11909 cells and 107829 ATAC peaks
```

``` r
atac_sc
#> <ScATAC> object with 11909 cells and 107829 ATAC peaks from hg38.
#> id: "PBMC"
#> description: "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#> loaded from: 'pbmc_data/matrix.mtx'
#> Slots include:
#>   • `@mat`: a numeric matrix where rows are peaks and columns are cells. Can be a sparse matrix.
#>   • `@peaks`: a misha intervals set with the peak definitions.
#>   • `@genome`: genome assembly of the peaks
```

### Filter peaks by coverage and/or length

``` r
atac_sc <- filter_features(scatac = atac_sc, minimal_max_umi = 3, min_peak_length = 200, max_peak_length = 1000)
```

### Project RNA metacells

``` r
data(cell_to_metacell_pbmc_example)
head(cell_to_metacell_pbmc_example)
#> # A tibble: 6 x 2
#>              cell_id metacell
#> 1 AAACAGCCAATCCCTT-1       44
#> 2 AAACAGCCAATGCGCT-1       22
#> 3 AAACAGCCACCAACCG-1        7
#> 4 AAACAGCCAGGATAAC-1       24
#> 5 AAACAGCCAGTTTACG-1       32
#> 6 AAACATGCAAGGTCCT-1       30
```

``` r
atac_mc <- project_atac_on_mc(atac_sc, cell_to_metacell_pbmc_example)
#> ℹ 3198 cells (out of 11909) do not have a metacell and have been removed.
#> ℹ Removed 142 all-zero peaks
#> • Setting egc cell size to 939452.6 (the 0.1 quantile of metacell sizes)
#> ✔ Created a new McATAC object with 97 metacells and 107687 ATAC peaks.
atac_mc
#> <McATAC> object with 97 metacells and 107687 ATAC peaks from hg38.
#> id: "PBMC"
#> description: "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#> Slots include:
#>   • `@mat`: a numeric matrix where rows are peaks and columns are metacells. Can be a sparse matrix.
#>   • `@peaks`: a misha intervals set with the peak definitions.
#>   • `@genome`: genome assembly of the peaks
#>   • `@egc`: a numeric matrix which contains normalized metacell accessibility.
#>   • `@fp`: a matrix showing for each peak (row) the relative enrichment of umis in log2 scale.
```

``` r
atac_mc
#> <McATAC> object with 97 metacells and 107687 ATAC peaks from hg38.
#> id: "PBMC"
#> description: "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#> Slots include:
#>   • `@mat`: a numeric matrix where rows are peaks and columns are metacells. Can be a sparse matrix.
#>   • `@peaks`: a misha intervals set with the peak definitions.
#>   • `@genome`: genome assembly of the peaks
#>   • `@egc`: a numeric matrix which contains normalized metacell accessibility.
#>   • `@fp`: a matrix showing for each peak (row) the relative enrichment of umis in log2 scale.
```

See more at the [vignette](https://tanaylab.github.io/mcATAC/articles/mcATAC.html)

### Add metadata

``` r
data(mcmd)
atac_mc <- add_mc_metadata(atac_mc, mcmd)
```

## Import RNA expression data

``` r
data(rna_mc_mat)
atac_mc <- add_mc_rna(atac_mc, rna_mc_mat)
```

``` r
plot_atac_rna(atac_mc, "CD4")
#> → The gene "CD4" has 9 alternative promoters. Summing the ATAC signal from all of them.
#> → The gene "CD4" has multiple (3) peaks within 500 bp of its TSS. Summing the ATAC signal from all of them.
```

<img src="man/figures/README-atac-rna-scatter-1-1.png" width="100%" />
