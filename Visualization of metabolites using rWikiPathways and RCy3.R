#Visualization of metabolite accumulation using pathway taken from Wikipathways
#Ken Kamiya, wrote at 2019/06/25

#To visualize the increase or decrease of the metabolites using pathway in WikiPathways.

#I used these assets.
#Pathway -> AtMetExpress Overview - Wikipathways https://www.wikipathways.org/index.php/Pathway:WP3622
#Metabolite data -> AraMetLeaves - attached to DiffCorr package
#The Arabidopsis metabolome of the aerial regions of individual WT plants and the mto1 and tt4 mutants were analyzed by GC-TOF/MS.
#mto1 shows strong methionine accumulation.
#DiffCorr -> https://cran.r-project.org/package=DiffCorr
#Reference of metabolite data -> M. Kusano et al, (2007) https://doi.org/10.1186/1752-0509-1-53


### Step 0. Install Cytoscape ###
#Before you run this script, you must install Cytoscape.
# https://cytoscape.org/download.html
#If you haven't, install this app from this link.


#### From here, please run with Cytoscape running ####


### Step 1. Prepare packages and dataset ###

#Prepare these packages and install "Wikipathways" app to Cytoscape
if(!"rWikiPathways" %in% installed.packages()){
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("rWikiPathways")
    BiocManager::install("RCy3")
}
library(rWikiPathways)
library(RCy3)

RCy3::installApp('WikiPathways')  #only available in Cytoscape 3.7.0 and above

#Prepare MetaboAnalystR
install.packages("pacman")
library(pacman)
pacman::p_load(Rserve, ellipse, scatterplot3d, Cairo, randomForest, caTools, e1071, som, impute, pcaMethods, RJSONIO, ROCR, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, pheatmap, SSPA, sva, Rcpp, pROC, data.table, limma, car, fitdistrplus, lars, Hmisc, magrittr, methods, xtable, pls, caret, lattice, igraph, gplots, KEGGgraph, reshape, RColorBrewer, tibble, siggenes, plotly, xcms, CAMERA, fgsea, MSnbase, BiocParallel, metap, reshape2, scales)
install.packages("devtools")
library(devtools)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))
library(MetaboAnalystR)

#Prepare metabolite data
install.packages("DiffCorr")
library(DiffCorr)
data("AraMetLeaves") #data for visualization



### Step 2. Gather information necessary for visualization ###

#Get graphIds from focused pathway
#https://www.wikipathways.org/index.php/Help:WikiPathways_Webservice/API#GraphId
pathway.ChEBIIDs <- rWikiPathways::getXrefList(pathway="WP3622", systemCode="Ce")
pathway.graphIds = NULL
for(i in 1:length(pathway.ChEBIIDs)){
  pathways <- rWikiPathways::findPathwaysByXref(pathway.ChEBIIDs[[i]], systemCode="Ce") 
  targeted.pathway <- pathways[purrr::map_lgl(pathways, ~ .$id == "WP3622")][[1]]
  pathway.graphIds[[i]] <- as.character(targeted.pathway$fields$graphId$values)
}


#Get metabolite ids lists using MetaboAnalystR
mSet<-InitDataObjects("NA", "utils", FALSE)
cmpd.vec<- rownames(AraMetLeaves)
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name", T, T, T, T, T);
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "Inositol-1-phosphate");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "Inositol-1-phosphate", "Inositol phosphate");
mSet<-PerformDetailMatch(mSet, "D-Glucose-6-phosphate");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "D-Glucose-6-phosphate", "Beta-D-Glucose 6-phosphate");

#Extract metabolites with graphID in the pathway
mSet.map.table = mSet$dataSet$map.table
mSet.graphIds = NULL
for(i in 1:length(mSet.map.table[,"ChEBI"])){
  pathways <- rWikiPathways::findPathwaysByXref(mSet.map.table[,"ChEBI"][[i]], systemCode="Ce")
  targeted.pathway = pathways[purrr::map_lgl(pathways, ~ .$id == "WP3622")]
  if(length(targeted.pathway)==0){
    mSet.graphIds[[i]] = "NA"
  }
  if(!length(targeted.pathway)==0){
    mSet.graphIds[[i]] <- targeted.pathway[[1]]$fields$graphId$values
  }
}

#Calculate log2FC from Col0.1 and mto1.1
Col0.index <- grep("Col0", colnames(AraMetLeaves), perl = TRUE)
mto1.index <- grep("mto1", colnames(AraMetLeaves), perl = TRUE)
Col0Mean <-  apply(AraMetLeaves[, Col0.index], 1, mean)
mto1Mean <-  apply(AraMetLeaves[, mto1.index], 1, mean)
log2FC <-  log2(mto1Mean / Col0Mean)
min.log2FC <-  min(log2FC, na.rm = TRUE)
max.log2FC <-  max(log2FC, na.rm = TRUE)
abs.log2FC <-  max(abs(min.log2FC), max.log2FC)
data.values <-  c(-abs.log2FC, 0, abs.log2FC)
##Make matrix
mSet.map.table.graphId.log2FC = data.frame(mSet.map.table, mSet.graphIds, log2FC)



### Step 3. Connect to Cytoscape and visalize ###

#Connect to Cytoscape
RCy3::cytoscapePing()
#Toss the pathway to Cytoscape
RCy3::commandsRun('wikipathways import-as-pathway id=WP3622') 
RCy3::toggleGraphicsDetails()
#load metabolite data and visualize
RCy3::loadTableData(mSet.map.table.graphId.log2FC, data.key.column = "mSet.graphIds", table.key.column = "GraphID")
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
RCy3::setNodeColorMapping("log2FC", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")

#emphasize methionine
RCy3::setNodeBorderColorBypass(node.names = "Methionine", new.colors = "#FF0000")
RCy3::setNodeBorderWidthBypass(node.names = "Methionine", new.sizes = 7)

table = RCy3::getTableColumns()
conv.table= table$name
names(conv.table) = table$GraphID

#Gray out metabolite nodes without data
RCy3::setNodeColorBypass(node.names = conv.table[setdiff(pathway.graphIds, mSet.graphIds)], new.colors = "#c0c0c0")
