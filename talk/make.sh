#!/bin/bash
Rscript -e "library(knitr);knit('talk.Rmd')"
pandoc -t beamer talk.md -o talk.pdf
open talk.pdf
