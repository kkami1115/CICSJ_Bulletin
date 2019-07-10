Visualization of metabolites using rWikiPathways and RCy3
================
ken kamiya
2019/07/07

# Abstract

This is the support page for “メタボロームデータ解析および解釈に資する可視化手法” in CICSJ
Bulletin July 2019.

This page describes the script in detail regarding the last chapter of
pathway visualization.

# Data

Datasets and assets I used are here.

  - Pathway AtMetExpress Overview - Wikipathways
    <https://www.wikipathways.org/index.php/Pathway:WP3622>

  - Metabolite data AraMetLeaves - attached to DiffCorr package

The Arabidopsis metabolome of the aerial regions of individual WT plants
and the mto1 and tt4 mutants were analyzed by GC-TOF/MS. mto1 shows
strong methionine accumulation.

  - DiffCorr: <https://cran.r-project.org/package=DiffCorr>
  - Reference of metabolite data: M. Kusano et al, (2007)
    <https://doi.org/10.1186/1752-0509-1-53>

# Step 0: Intstall Cytoscape

If you haven’t yet installed Cytoscape, install it before running this
script.

<https://cytoscape.org/download.html>

The following script only works if you have launched Cytoscape.

# Step 1: Prepare R Packages and Metabolite Data

This pattern of installation had succeeded on my environment(Ubuntu
16.04, R 3.6.0).

It doesn’t matter how you do it, so install “rWikiPathways”, “RCy3”,
“MetabAnalystR” according to your own environment.

  - Especially, Installing MetaboAnalystR is extremely difficult, so
    please refer to (<https://github.com/xia-lab/MetabanaristR>) for
    installation.

<!-- end list -->

``` r
install.packages("BiocManager", repos = "https://cran.rstudio.com")
```

    ## Installing package into '/home/kamiya/R/x86_64-pc-linux-gnu-library/3.6'
    ## (as 'lib' is unspecified)

``` r
BiocManager::install("rWikiPathways")
```

    ## Bioconductor version 3.9 (BiocManager 1.30.4), R 3.6.1 (2019-07-05)

    ## Installing package(s) 'rWikiPathways'

``` r
library(rWikiPathways)
BiocManager::install("RCy3")
```

    ## Bioconductor version 3.9 (BiocManager 1.30.4), R 3.6.1 (2019-07-05)

    ## Installing package(s) 'RCy3'

    ## Warning in install.packages(pkgs = doing, lib = lib, repos = repos, ...):
    ## installation of package 'RCy3' had non-zero exit status

``` r
library(RCy3)
```

    ## Registered S3 method overwritten by 'R.oo':
    ##   method        from       
    ##   throw.default R.methodsS3

``` r
#Prepare MetaboAnalystR
install.packages("pacman", repos = "https://cran.rstudio.com")
```

    ## Installing package into '/home/kamiya/R/x86_64-pc-linux-gnu-library/3.6'
    ## (as 'lib' is unspecified)

``` r
library(pacman)
pacman::p_load(Rserve, ellipse, scatterplot3d, Cairo, randomForest, caTools, e1071, som, impute, pcaMethods, RJSONIO, ROCR, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, pheatmap, SSPA, sva, Rcpp, pROC, data.table, limma, car, fitdistrplus, lars, Hmisc, magrittr, methods, xtable, pls, caret, lattice, igraph, gplots, KEGGgraph, reshape, RColorBrewer, tibble, siggenes, plotly, xcms, CAMERA, fgsea, MSnbase, BiocParallel, metap, reshape2, scales)
install.packages("devtools", repos = "https://cran.rstudio.com")
```

    ## Installing package into '/home/kamiya/R/x86_64-pc-linux-gnu-library/3.6'
    ## (as 'lib' is unspecified)

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))
```

    ## Skipping install of 'MetaboAnalystR' from a github remote, the SHA1 (593819a8) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library(MetaboAnalystR)
```

    ## Warning: replacing previous import 'CAMERA::annotate' by
    ## 'ggplot2::annotate' when loading 'MetaboAnalystR'

    ## Warning: replacing previous import 'MSnbase::normalize' by
    ## 'igraph::normalize' when loading 'MetaboAnalystR'

    ## Warning: replacing previous import 'reshape::recast' by 'reshape2::recast'
    ## when loading 'MetaboAnalystR'

    ## Warning: replacing previous import 'reshape::melt' by 'reshape2::melt' when
    ## loading 'MetaboAnalystR'

    ## Warning: replacing previous import 'reshape::colsplit' by
    ## 'reshape2::colsplit' when loading 'MetaboAnalystR'

    ## Warning: replacing previous import 'igraph::normalize' by 'xcms::normalize'
    ## when loading 'MetaboAnalystR'

    ## Warning: replacing previous import 'igraph::groups' by 'xcms::groups' when
    ## loading 'MetaboAnalystR'

    ## 
    ## Attaching package: 'MetaboAnalystR'

    ## The following object is masked from 'package:caret':
    ## 
    ##     splsda

``` r
install.packages("DiffCorr", repos = "https://cran.rstudio.com")
```

    ## Installing package into '/home/kamiya/R/x86_64-pc-linux-gnu-library/3.6'
    ## (as 'lib' is unspecified)

``` r
library(DiffCorr)
```

    ## Loading required package: fdrtool

``` r
data("AraMetLeaves") #data for visualization

Col0.index <- grep("Col0", colnames(AraMetLeaves), perl = TRUE)
mto1.index <- grep("mto1", colnames(AraMetLeaves), perl = TRUE)
Col0Mean <-  apply(AraMetLeaves[, Col0.index], 1, mean)
mto1Mean <-  apply(AraMetLeaves[, mto1.index], 1, mean)
```

# Step 2: Gather informations (espescially metabolite ID)

## Get graphIds from the pathway

Before getting pathways you need to know “graphId” in
WikiPathways.

<https://www.wikipathways.org/index.php/Help:WikiPathways_Webservice/API#GraphId>

graphId is an individual ID of a metabolite assigned in the pathway. By
using this, the pathway metabolite and metabolite data are
linked.

``` r
pathway.HMDBIDs <- rWikiPathways::getXrefList(pathway="WP3622", systemCode="Ch")
pathway.graphIds <-  NULL
for(i in 1:length(pathway.HMDBIDs)){
  pathways <- rWikiPathways::findPathwaysByXref(pathway.HMDBIDs[[i]], systemCode="Ch") 
  targeted.pathway <- pathways[purrr::map_lgl(pathways, ~ .$id == "WP3622")][[1]]
  pathway.graphIds[[i]] <- as.character(targeted.pathway$fields$graphId$values)
}
```

## Get IDs of metabolite data using MetaboAnalystR

Metabolite data names are written in various ways and are not stable.

So we convert names to IDs such as HMDB or KEGG using MetaboAnalystR.

The same can be done with MetaboAnalyst (<https://metaboanalyst.ca>) on
the Web.

Entering a metabolite name returns the match result recognized by
MetaBoAnalyst, so you can absorb minor differences in metabolite names.

In addition, MetaboAnalyst can analyze while generating R code, so it is
also possible to use the displayed R code in MetaBoAnalystR after
performing ID conversion manually with MetaboAnalyst on the Web.

``` r
mSet<-InitDataObjects("NA", "utils", FALSE)
```

    ## [1] "MetaboAnalyst R objects initialized ..."

``` r
cmpd.vec<- rownames(AraMetLeaves)
mSet<-Setup.MapData(mSet, cmpd.vec)
mSet<-CrossReferencing(mSet, "name", T, T, T, T, T)
```

    ## [1] "Loaded files from MetaboAnalyst web-server."
    ## [1] "Loaded files from MetaboAnalyst web-server."

``` r
mSet<-CreateMappingResultTable(mSet)
```

    ## [1] "Loaded files from MetaboAnalyst web-server."

``` r
mSet<-PerformDetailMatch(mSet, "Inositol-1-phosphate")
```

    ## [1] "Loaded files from MetaboAnalyst web-server."
    ## [1] "Loaded files from MetaboAnalyst web-server."

``` r
mSet<-GetCandidateList(mSet)
```

    ## [1] "Loaded files from MetaboAnalyst web-server."

``` r
mSet<-SetCandidate(mSet, "Inositol-1-phosphate", "Inositol phosphate")
```

    ## [1] "Loaded files from MetaboAnalyst web-server."

``` r
mSet<-PerformDetailMatch(mSet, "D-Glucose-6-phosphate")
```

    ## [1] "Loaded files from MetaboAnalyst web-server."
    ## [1] "Loaded files from MetaboAnalyst web-server."

``` r
mSet<-GetCandidateList(mSet)
```

    ## [1] "Loaded files from MetaboAnalyst web-server."

``` r
mSet<-SetCandidate(mSet, "D-Glucose-6-phosphate", "Beta-D-Glucose 6-phosphate")
```

    ## [1] "Loaded files from MetaboAnalyst web-server."

## Get metabolites with graphId in the pathway

The metabolites having pathwayId in the metabolite data are extracted.

``` r
mSet.map.table <- mSet$dataSet$map.table
mSet.graphIds <- NULL
for(i in 1:length(mSet.map.table[,"HMDB"])){
  pathways <- findPathwaysByXref(mSet.map.table[,"HMDB"][[i]], systemCode="Ch")
  targeted.pathway <- pathways[purrr::map_lgl(pathways, ~ .$id == "WP3622")]
  if(length(targeted.pathway) == 0){
    mSet.graphIds[[i]] <- "NA"
  }
  if(!length(targeted.pathway) == 0){
    mSet.graphIds[[i]] <- targeted.pathway[[1]]$fields$graphId$values
  }
}
```

## Calculate log2FC from Col0 and mto1

For visualization, log2FC is calculated from the mean of Col0 and the
mean of mto1 mutants.

``` r
log2FC <-  log2(mto1Mean / Col0Mean)
min.log2FC <-  min(log2FC, na.rm = TRUE)
max.log2FC <-  max(log2FC, na.rm = TRUE)
abs.log2FC <-  max(abs(min.log2FC), max.log2FC)
data.values <-  c(-abs.log2FC, 0, abs.log2FC)
```

## Make final matrix

Create a matrix to pass to
Cytoscape.

``` r
mSet.map.table.graphId.log2FC <-  data.frame(mSet.map.table, mSet.graphIds, log2FC)
```

# Step 3: Connect to Cytoscape and Visualize

Connect to Cytoscape using RCy3 for visualization.

**Start Cytoscape before running this chapter**

``` r
#Install Wikipathways addin to Cytoscape
installApp('WikiPathways')  #only available in Cytoscape 3.7.0 and above
```

    ## [1] "App WikiPathways installed"

``` r
#Connect to Cytoscape
cytoscapePing()
```

    ## [1] "You are connected to Cytoscape!"

``` r
#Pass the pathway to Cytoscape
RCy3::commandsRun('wikipathways import-as-pathway id=WP3622') 
toggleGraphicsDetails()

#load metabolite data and visualize 
loadTableData(mSet.map.table.graphId.log2FC, data.key.column = "mSet.graphIds", table.key.column = "GraphID")
```

    ## [1] "Success: Data loaded in defaultnode table"

``` r
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2FC", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")

#emphasize methionine
setNodeBorderColorBypass(node.names = "Methionine", new.colors = "#FF0000")
setNodeBorderWidthBypass(node.names = "Methionine", new.sizes = 10)
```
