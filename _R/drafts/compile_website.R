#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================

library(rmarkdown)
library(here)

# Compile site
render_site(here())

# Compile pdf
# pdf
# options(knitr.table.format="latex")
# pdf_format <- bookdown::pdf_document2(citation_package="natbib", fig_caption=TRUE, keep_tex=FALSE, keep_md=FALSE,
# 																			latex_engine="pdflatex", number_sections=TRUE, toc=FALSE)
# bookdown::render_book(here("index.Rmd"), output_format=pdf_format)

# EOF