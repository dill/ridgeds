#!/bin/bash
Rscript -e "library(knitr);knit('talk.Rmd')"
pandoc -t beamer talk.md -H head.tex -o talk.pdf
rm talk.md
open talk.pdf
