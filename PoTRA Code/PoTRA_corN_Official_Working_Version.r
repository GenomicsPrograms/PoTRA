

PoTRA.corN <- function(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database, PR.quantile) {
  
  require(BiocGenerics)
  require(graph)
  require(graphite)
  require(igraph)  
  require(org.Hs.eg.db) 
  Fishertest<-c()
  TheNumOfHubGene.normal<-c()
  TheNumOfHubGene.case<-c()
  E.normal<-c()
  E.case<-c()
  length.pathway<-c()
  kstest<-c()
  pathwaynames <- c()
 
  humanReactome <- pathways("hsapiens", "reactome")
  humanBiocarta <- pathways("hsapiens", "biocarta")
  humanKEGG <- pathways("hsapiens", "kegg")
  
   for (x in 1:length(Pathway.database[1:length(humanKEGG)])){
   
    print(x)
    p0 <-Pathway.database[[x]]
    pathwaynames[x] <- p0@title
    p <- convertIdentifiers(p0, "entrez")
    g<-pathwayGraph(p) 
    nodelist<-nodes(g)
    graph.path<-igraph.from.graphNEL(g)    
	## Plot graph
	#plot(graph.path, edge.arrow.size=.5, vertex.color="gold", vertex.size=5, 
    #vertex.frame.color="gray", vertex.label.color="black", 
    #vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.2) 
	
	
	genelist_reformatted <- sprintf('ENTREZID:%s', genelist) 
	length.intersect<-length(intersect(unlist(nodelist),unlist(genelist_reformatted)))
    length.pathway[x]<-length.intersect
    graph.path<-induced_subgraph(graph.path, as.character(intersect(unlist(nodelist),unlist(genelist_reformatted))))
    #plot(graph.path)	
	
	
    if (length.intersect<5){
      next
    }else{
      
      #collect expression data of genes for a specific pathway across normal and tumor samples.
      
      path<-data.frame(matrix(0,length.intersect,(Num.sample.normal+Num.sample.case)))
      a<- c()
      
      for (j in 1:length.intersect){
        a[j]<-intersect(unlist(nodelist),unlist(genelist_reformatted))[j]  
        
        path[j,]<-mydata[which(genelist_reformatted==a[j]),]  #collect expression data of genes for a specific pathway across normal and tumor samples.
      }
      
      ##Construct a gene-gene network for normal samples and calculate PageRank values for each gene in this network.
      
      cor.normal <- apply(path[,1:Num.sample.normal], 1, function(x) { apply(path[,1:Num.sample.normal], 1, function(y) { cor.test(x,y)[[3]] })})
      
      cor.normal<-as.matrix(cor.normal) 
      
      cor.normal.adj<-matrix(p.adjust(cor.normal,method="fdr"),length.intersect,length.intersect)
      
      cor.normal.adj[ cor.normal.adj > 0.05 ] <- 0
      cor.normal.adj[ is.na(cor.normal.adj)] <- 0
      diag(cor.normal.adj) <- 0
      graph.normal<-graph.adjacency(cor.normal.adj,weighted=TRUE,mode="undirected")
      E.normal[x]<-length(E(graph.normal))
      
      PR.normal<-page.rank(graph.normal,direct=FALSE)$vector
      
      ##Construct a gene-gene network for tumor samples and calculate PageRank values for each gene in this network.
      
      cor.case <- apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(x) { apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(y) { cor.test(x,y)[[3]] })})
      
      cor.case<-as.matrix(cor.case) 
      cor.case.adj<-matrix(p.adjust(cor.case,method="fdr"),length.intersect,length.intersect)
      
      cor.case.adj[ cor.case.adj > 0.05 ] <- 0
      cor.case.adj[ is.na(cor.case.adj)] <- 0
      diag(cor.case.adj) <- 0
      graph.case<-graph.adjacency(cor.case.adj,weighted=TRUE,mode="undirected")
      E.case[x]<-length(E(graph.case))
      
      PR.case<-page.rank(graph.case,direct=FALSE)$vector
	 
 
      matrix.HCC<- matrix("NA",2*length.intersect,2)
      rownames(matrix.HCC)<-as.character(c(PR.normal,PR.case))
      colnames(matrix.HCC)<-c("Disease_status","PageRank")
      
      matrix.HCC[,1]<-c(rep("Normal",length.intersect), rep("Cancer",length.intersect))
      
      loc.largePR<-which(as.numeric(rownames(matrix.HCC))>=quantile(PR.normal,PR.quantile))
      loc.smallPR<-which(as.numeric(rownames(matrix.HCC))<quantile(PR.normal,PR.quantile))
      
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
	  tryCatch(
      kstest[x]<-ks.test(PR.normal,PR.case)$p.value, 
	  error=function(e) print("."))  
	
      
      if (E.normal[x]<E.case[x]){
        Fishertest[x]<-1
        kstest[x]<-1
      }else{
        Fishertest[x]<-Fishertest[x]
        kstest[x]<-kstest[x]
      }
       
    }
  }
  return(list(Fishertest.p.value=Fishertest,KStest.p.value=kstest,LengthOfPathway=length.pathway,TheNumberOfHubGenes.normal=TheNumOfHubGene.normal,TheNumOfHubGene.case=TheNumOfHubGene.case,TheNumberOfEdges.normal=E.normal,TheNumberOfEdges.case=E.case,PathwayName=pathwaynames))
}



mydata <- read.table(args[1], sep="\t", header=TRUE)
genes <- as.data.frame(rownames(mydata))
names(genes) <- c("entrez")
genelist <- genes[,1]

library(org.Hs.eg.db) 
library(BiocGenerics) 
library(graph) 
library(graphite) 
library(igraph) 

humanKEGG <- pathways("hsapiens", "kegg")
  
Pathway.database <- humanKEGG
PR.quantile = as.numeric(args[2])
Num.sample.normal= as.numeric(args[3]) ## update value, value must be at least 20
Num.sample.case=as.numeric(args[4])  ## update value, value must be at least 20

results.cor <-PoTRA.corN(mydata=mydata,genelist=genelist,Num.sample.normal=Num.sample.normal,Num.sample.case=Num.sample.case,Pathway.database=Pathway.database,PR.quantile=PR.quantile)



#######################Post-processing########################################


fishertest <- as.data.frame(unlist(results.cor[1]))
names(fishertest) <- c("fishertest.pvalue")
rownames(fishertest) <- NULL

kstest <- as.data.frame(unlist(results.cor[2]))
names(kstest) <- c("kstest.pvalue")
rownames(kstest) <- NULL

length.pathway <- as.data.frame(unlist(results.cor[3]))
names(length.pathway) <- c("length.pathway")
rownames(length.pathway) <- NULL


TheNumberOfHubGenes.normal <- as.data.frame(unlist(results.cor[4]))
names(TheNumberOfHubGenes.normal) <- c("TheNumberOfHubGenes.normal")
rownames(TheNumberOfHubGenes.normal) <- NULL


TheNumberOfHubGenes.case <- as.data.frame(unlist(results.cor[5]))
names(TheNumberOfHubGenes.case) <- c("TheNumberOfHubGenes.case")
rownames(TheNumberOfHubGenes.case) <- NULL


TheNumberOfEdges.normal <- as.data.frame(unlist(results.cor[6]))
names(TheNumberOfEdges.normal) <- c("TheNumberOfEdges.normal")
rownames(TheNumberOfEdges.normal) <- NULL


TheNumberOfEdges.case <- as.data.frame(unlist(results.cor[7]))
names(TheNumberOfEdges.case) <- c("TheNumberOfEdges.case")
rownames(TheNumberOfEdges.case) <- NULL

Pathway <- as.data.frame(unlist(results.cor[8]))
names(Pathway) <- c("Pathway")
rownames(Pathway) <- NULL

## Remove trailing rows in Pathway 

diff = -1*(nrow(Pathway) - nrow(fishertest))
Pathways = head(Pathway, diff)
Pathway.Length = head(length.pathway, diff)

## Sort by fishertest ascending and then for the Pathway var, apply an alpha sort
df_result.corN <- data.frame(Pathways, fishertest, kstest, TheNumberOfHubGenes.normal, TheNumberOfHubGenes.case, TheNumberOfEdges.normal, TheNumberOfEdges.case, Pathway.Length)


results0 <- df_result.corN[order(fishertest),] 
results0$Rank <- 1:nrow(results0) 

potra_results <- results0[order(results0$Pathway),]

data1 <- data.frame(lapply(potra_results, function(x) { gsub(", ", "_", x)}))
data2 <- data.frame(lapply(data1, function(x) { gsub(" ", "_", x)}))

write.table(data2, "PoTRA_Results.txt", sep="\t", quote=FALSE, row.name=FALSE)

