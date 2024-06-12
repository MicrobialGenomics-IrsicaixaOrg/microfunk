
<!-- README.md is generated from README.Rmd. Please edit that file -->

# microfunk

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/MicrobialGenomics-IrsicaixaOrg/microfunk)](https://github.com/MicrobialGenomics-IrsicaixaOrg/microfunk/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/MicrobialGenomics-IrsicaixaOrg/microfunk)](https://github.com/MicrobialGenomics-IrsicaixaOrg/microfunk/pulls)
<!-- badges: end -->

The goal of `microfunk` is to …

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `microfunk` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("microfunk")
```

And the development version from
[GitHub](https://github.com/juditfarreb/microfunk) with:

``` r
BiocManager::install("juditfarreb/microfunk")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("microfunk")
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!

## Citation

Below is the citation output from using `citation('microfunk')` in R.
Please run this yourself to check for any updates on how to cite
**microfunk**.

``` r
print(citation('microfunk'), bibtex = TRUE)
#> Warning in citation("microfunk"): could not determine year for 'microfunk' from
#> package DESCRIPTION file
#> To cite package 'microfunk' in publications use:
#> 
#>   Farré-Badia J (????). _microfunk: What the Package Does (One Line,
#>   Title Case)_. R package version 0.99.0,
#>   https://microbialgenomics-irsicaixaorg.github.io/microfunk/,
#>   <https://github.com/MicrobialGenomics-IrsicaixaOrg/microfunk>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {microfunk: What the Package Does (One Line, Title Case)},
#>     author = {Judit Farré-Badia},
#>     note = {R package version 0.99.0, 
#> https://microbialgenomics-irsicaixaorg.github.io/microfunk/},
#>     url = {https://github.com/MicrobialGenomics-IrsicaixaOrg/microfunk},
#>   }
```

Please note that the `microfunk` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `microfunk` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools
