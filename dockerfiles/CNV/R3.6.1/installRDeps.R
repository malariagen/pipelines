#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos=r)
})


## Install Biostrings
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")

## Old Code
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("BiocParallel")

library(BiocManager)
BiocManager::valid()

print("checking Biostrings")

if ( ! library("Biostrings", character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }


## Install other libraries
install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
install.packages("stringr", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
# install.packages("parallel", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

print("checking stringi")
if ( ! library('stringi', character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }
print("checking stringr")
if ( ! library("stringr", character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }
print("checking parallel")
if ( ! library("parallel", character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }