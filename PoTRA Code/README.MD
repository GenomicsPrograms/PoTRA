
# Usage <br>

The R script to perform the main part (Fisher’s exact test) of PoTRA, the alternative part (K.S. test) of PoTRA, and the overlay plot of kernel density for the PR scores: <b>PoTRA.R</b><br>

Example usage of PoTRA.R: <b>PoTRA-example.R</b><br>

Data used in PoTRA-example.R: <b>PoTRA-example-data.Rdata</b><br>


library(BiocGenerics) <br>
library(graph) <br>
library(graphite) <br>
library(igraph) <br>

<b>Choose your database before running the PoTRA.cor commandline: <br> </b>
<code>humanKEGG <- pathways("hsapiens", "kegg") </code>


<b>For PoTRA (Fisher’s exact test and K.S. test) based on correlation networks:</b> <br>

<code>PoTRA.cor <- function(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database, PR.quantile)</code>
 
<b>Example Usage: </b><br>
<code>results.cor <-PoTRA.cor(mydata=mydata,genelist=genelist,Num.sample.normal=8,Num.sample.case=8,Pathway.database=humanKEGG[1:10],PR.quantile=0.95)</code>


<b>For PoTRA (Fisher’s exact test and K.S. test) based on combined networks:</b> <br>

<code>PoTRA.comb <- function(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database)</code>

<b>Example Usage: </b><br>
<code>results.comb <-PoTRA.comb(mydata=mydata,genelist=genelist,Num.sample.normal=8,Num.sample.case=8,Pathway.database=humanKEGG[1:10],PR.quantile=0.95)</code>


<b>For the kernel density plot:</b> <br>

<code>results.plot <- overlayplot(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database)</code>

 
<b>Example Usage: </b><br>

<code>results.plot <- overlayplot(mydata=mydata,genelist=genelist,Num.sample.normal=8,Num.sample.case=8,Pathway.database=humanKEGG[[1]])</code>

# Arguments

<html>
<body>

<table>
   <tr>
    <td>mydata</td>
    <td>A gene expression dataset (matrix). Rows represent genes, and columns represent samples (from control to case)</td>
   </tr>
   <tr>
    <td>genelist</td>
    <td>A list of gene names</td>
   </tr>
    <tr>
    <td>Num.sample.normal</td>
    <td>The number of normal samples</td>
   </tr>
   <tr>
    <td>Pathway.database</td>
    <td>The pathway database used in the study, such as, KEGG, Biocarta and Reactome</td>
   </tr>
    <tr>
    <td>PR.quantile</td>
    <td>The percentile of PageRank scores as a cutoff for hub genes </td>
    </tr>
</table>

</body>
</html>

   
# Values

 
<html>
<body>

<table>
   <tr>
    <td>Fishertest.p.value</td>
    <td>The p-value of the Fisher’s exact test</td>
   </tr>
   <tr>
    <td>KStest.p.value</td>
    <td>The p-value of the K.S. test</td>
   </tr>
    <tr>
    <td>LengthOfPathway</td>
    <td>The length of pathways</td>
   </tr>
   <tr>
    <td>TheNumberOfHubGenes.normal</td>
    <td>the number of hub genes for normal samples</td>
   </tr>
    <tr>
    <td>PR.quantile</td>
    <td>The percentile of PageRank scores as a cutoff for hub genes </td>
    </tr>
    <tr>
    <td>TheNumOfHubGene.case</td>
    <td>The number of hub genes for cancer samples </td>
    </tr>
    <tr>
    <td>TheNumberOfEdges.normal</td>
    <td>The number of edges in the network for normal samples</td>
    </tr>
    <tr>
    <td>TheNumberOfEdges.case</td>
    <td>The number of edges in the network for cancer samples </td>
    </tr>
</table>

</body>
</html>
 


<b>Examples </b>

See PoTRA-example.R

 
# Credits

Authors: Chaoxing Li, Li Liu and Valentin Dinu

 

<b>If you use or modify the code, please cite:</b>

“Pathways of Topological Rank Analysis (PoTRA): A Novel Method to Detect Pathways Involved in Cancer”, 2018.

 

<b>Issues?</b>

Please email Chaoxing Li <chaoxing@asu.edu> or Valentin Dinu <Valentin.Dinu@asu.edu> if you have any questions, concerns or feedback.
