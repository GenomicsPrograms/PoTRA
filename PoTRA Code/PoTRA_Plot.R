#!/usr/bin/env Rscript
## ML: Added code so that the program can be called via the commandline
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please provide at least one argument (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[6] = "Results_Overlay_Plot.png"
}
## Usage: Rscript mydata genelist Num.sample.normal Num.sample.case Pathway.database Results_Overlay_Plot.png

## args[1] = mydata
## args[2] = genelist
## args[3] = Num.sample.normal
## args[4] = Num.sample.case
## args[5] = Pathway.database

overlayplot<-function(args[1],args[2],args[3],args[4],args[5]) {
  
  require(BiocGenerics)
  require(graph)
  require(graphite)
  require(igraph)
  
  
  mydata <- args[1]
  genelist <- args[2]
  Num.sample.normal <- args[3]
  Num.sample.case <- args[4]
  Pathway.database <- args[5]
  
  length.pathway<-c()
  
  humanReactome <- pathways("hsapiens", "reactome")
  humanBiocarta <- pathways("hsapiens", "biocarta")
  humanKEGG <- pathways("hsapiens", "kegg")
  
  plot.multi.dens <- function(s)
  {
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s))
    {
      junk.x = c(junk.x, density(s[[i]])$x)
      junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = "")
    for(i in 1:length(s))
    {
      lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
    }
  }
  
  p0 <-Pathway.database
  p <- convertIdentifiers(p0, "entrez")
  g <- pathwayGraph(p) 
  nodelist <- nodes(g)
  graph.path <- igraph.from.graphNEL(g)    
  length.intersect <- length(intersect(unlist(nodelist),unlist(genelist)))
  
  length.pathway<-length(nodelist)
  
  #collect expression data of genes for a specific pathway across normal and tumor samples.
  
  path<-data.frame(matrix(200,length.intersect,(Num.sample.normal+Num.sample.case)))
  a<- c()
  
  for (j in 1:length.intersect){
    a[j]<-intersect(unlist(nodelist),unlist(genelist))[j]  
    
    path[j,]<-mydata[which(genelist==a[j]),]  #collect expression data of genes for a specific pathway across normal and tumor samples.
  }
  
  ##Construct a gene-gene network for normal samples and calculate PageRank values for each gene in this network.
  
  cor.normal <- apply(path[,1:Num.sample.normal], 1, function(x) { apply(path[,1:Num.sample.normal], 1, function(y) { cor.test(x,y)[[3]] })})
  
  cor.normal<-as.matrix(cor.normal) 
  
  cor.normal.adj<-matrix(p.adjust(cor.normal,method="fdr"),length.intersect,length.intersect)
  
  cor.normal.adj[ cor.normal.adj > 0.05 ] <- 0
  cor.normal.adj[ is.na(cor.normal.adj)] <- 0
  diag(cor.normal.adj) <- 0
  graph.normal<-graph.adjacency(cor.normal.adj,weighted=TRUE,mode="undirected")
  
  PR.normal<-page.rank(graph.normal,direct=FALSE)
  
  ##Construct a gene-gene network for tumor samples and calculate PageRank values for each gene in this network.
  
  cor.case <- apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(x) { apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(y) { cor.test(x,y)[[3]] })})
  
  cor.case<-as.matrix(cor.case) 
  cor.case.adj<-matrix(p.adjust(cor.case,method="fdr"),length.intersect,length.intersect)
  
  cor.case.adj[ cor.case.adj > 0.05 ] <- 0
  cor.case.adj[ is.na(cor.case.adj)] <- 0
  diag(cor.case.adj) <- 0
  graph.case<-graph.adjacency(cor.case.adj,weighted=TRUE,mode="undirected")
  
  PR.case<-page.rank(graph.case,direct=FALSE)
  
  PoTRA.plot<-plot.multi.dens(list(PR.normal$vector,PR.case$vector))
}
