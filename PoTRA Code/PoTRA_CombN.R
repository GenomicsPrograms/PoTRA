#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please provide at least one argument (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[7] = "Results_CombN.txt"
}
## Usage: Rscript mydata genelist Num.sample.normal Num.sample.case Pathway.database PR.quantile Results_CombN.txt

## args[1] = mydata
## args[2] = genelist
## args[3] = Num.sample.normal
## args[4] = Num.sample.case
## args[5] = Pathway.database
## args[6] = PR.quantile

PoTRA.combN <- function(args[1], args[2], args[3], args[4], args[5], args[6]) {
  
  require(BiocGenerics)
  require(graph)
  require(graphite)
  require(igraph)
  
  
  mydata <- args[1]
  genelist <- args[2]
  Num.sample.normal <- args[3]
  Num.sample.case <- args[4]
  Pathway.database <- args[5]
  PR.quantile <- args[6]
  
  
  Fishertest<-c()
  TheNumOfHubGene.normal<-c()
  TheNumOfHubGene.case<-c()
  E.normal<-c()
  E.case<-c()
  length.pathway<-c()
  E.union.normal<-c()
  E.union.case<-c()
  kstest<-c()
  pathwaynames <- c()
  
  humanReactome <- pathways("hsapiens", "reactome")
  humanBiocarta <- pathways("hsapiens", "biocarta")
  humanKEGG <- pathways("hsapiens", "kegg")
  
  for (x in 1:length(Pathway.database)){
    print(x)
    p0 <-Pathway.database[[x]]
    pathwaynames[x] <- p0@title
    p <- convertIdentifiers(p0, "entrez")
    g<-pathwayGraph(p) 
    nodelist<-nodes(g)
    graph.path<-igraph.from.graphNEL(g)
    graph.path<-as.undirected(graph.path)    
    length.intersect<-length(intersect(unlist(nodelist),unlist(genelist)))
    
    length.pathway[x]<-length.intersect
    
    graph.path<-induced_subgraph(graph.path, as.character(intersect(unlist(nodelist),unlist(genelist))))
    
    if (length.intersect<5){
      next
    }else{
      
      #collect expression data of genes for a specific pathway across normal and tumor samples.
      
      path<-data.frame(matrix(0,length.intersect,(Num.sample.normal+Num.sample.case)))
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
      
      colnames(cor.normal.adj)<-intersect(unlist(nodelist),unlist(genelist))
      rownames(cor.normal.adj)<-intersect(unlist(nodelist),unlist(genelist))
      
      graph.normal<-graph.adjacency(cor.normal.adj,weighted=TRUE,mode="undirected")
      E.normal[x]<-length(E(graph.normal))
      PR.normal<-page.rank(graph.normal,direct=FALSE)$vector
      
      graph.union.normal<- intersection(graph.normal,graph.path)
      E.union.normal[x]<- length(E(graph.union.normal))
      PR.union.normal<-page.rank(graph.union.normal,direct=FALSE)$vector
      
      ##Construct a gene-gene network for tumor samples and calculate PageRank values for each gene in this network.
      
      cor.case <- apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(x) { apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(y) { cor.test(x,y)[[3]] })})
      
      cor.case<-as.matrix(cor.case) 
      cor.case.adj<-matrix(p.adjust(cor.case,method="fdr"),length.intersect,length.intersect)
      
      cor.case.adj[ cor.case.adj > 0.05 ] <- 0
      cor.case.adj[ is.na(cor.case.adj)] <- 0
      diag(cor.case.adj) <- 0
      
      colnames(cor.case.adj)<-intersect(unlist(nodelist),unlist(genelist))
      rownames(cor.case.adj)<-intersect(unlist(nodelist),unlist(genelist))
      
      graph.case<-graph.adjacency(cor.case.adj,weighted=TRUE,mode="undirected")
      E.case[x]<-length(E(graph.case))
      PR.case<-page.rank(graph.case,direct=FALSE)$vector
      
      graph.union.case<- intersection(graph.case,graph.path)
      E.union.case[x]<- length(E(graph.union.case))
      PR.union.case<-page.rank(graph.union.case,direct=FALSE)$vector
      ############################
      matrix.HCC<- matrix("NA",2*length.intersect,2)
      rownames(matrix.HCC)<-as.character(c(PR.union.normal,PR.union.case))
      colnames(matrix.HCC)<-c("Disease_status","PageRank")
      
      matrix.HCC[,1]<-c(rep("Normal",length.intersect), rep("Cancer",length.intersect))
      
      loc.largePR<-which(as.numeric(rownames(matrix.HCC))>=quantile(PR.union.normal,PR.quantile))
      loc.smallPR<-which(as.numeric(rownames(matrix.HCC))<quantile(PR.union.normal,PR.quantile))
      
      matrix.HCC[loc.largePR,2]<-"large_PageRank"
      matrix.HCC[loc.smallPR,2]<-"small_PageRank"
      
      table.HCC<-list(1,2)
      names(table.HCC)<-c("Disease_status","PageRank")
      
      table.HCC$Disease_status<-matrix("NA",2*length.intersect,2)
      table.HCC$PageRank<-matrix("NA",2*length.intersect,2)
      
      table.HCC$Disease_status<-matrix.HCC[,1]
      table.HCC$PageRank<-matrix.HCC[,2]
      
      cont.HCC<-table(table.HCC$Disease_status,table.HCC$PageRank)
      TheNumOfHubGene.normal[x]<-cont.HCC[2]
      TheNumOfHubGene.case[x]<-cont.HCC[1]
      
      if (dim(cont.HCC)[1]!=dim(cont.HCC)[2]){
        Fishertest[x]<-1
      }else{
        Fishertest[x]<-fisher.test(cont.HCC)$p.value
      }
      
      kstest[x]<-ks.test(PR.union.normal,PR.union.case)$p.value
      
      if (E.union.normal[x]<E.union.case[x]){
        Fishertest[x]<-1
        kstest[x]<-1
      }else{
        Fishertest[x]<-Fishertest[x]
        kstest[x]<-kstest[x]
      }
      ############################################
    }
  }
  return(list(Fishertest.p.value=Fishertest,KStest.p.value=kstest,LengthOfPathway=length.pathway,TheNumberOfHubGenes.normal=TheNumOfHubGene.normal,TheNumOfHubGene.case=TheNumOfHubGene.case,TheNumberOfEdges.normal=E.union.normal,TheNumberOfEdges.case=E.union.case,PathwayName=pathwaynames))
}
