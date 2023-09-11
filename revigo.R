library(ggplot2)
library(scales)
library(rrvgo)

go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")


heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

## -----------------------------------------------------------------------------
scatterPlot(simMatrix, reducedTerms)

## ---- eval=FALSE--------------------------------------------------------------
treemapPlot(reducedTerms)

## -----------------------------------------------------------------------------
wordcloudPlot(reducedTerms, min.freq=1, colors="black")

## ---- eval=FALSE--------------------------------------------------------------
#  rrvgo::shiny_rrvgo()

## ---- eval=FALSE--------------------------------------------------------------
#  my_new_fancy_orgdb_object <- 'org.Zz.eg.db'
#  hsGO <- GOSemSim::godata(my_new_fancy_orgdb_object, ont="MF")

## ----citation-----------------------------------------------------------------
citation("rrvgo")

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()
