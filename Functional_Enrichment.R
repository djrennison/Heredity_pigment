library(topGO)
library(Rgraphviz)

geneID2GO <- readMappings(file = "stickleGO.txt") 
geneUniverse <- names(geneID2GO)


#Chromosome 8 
candidate_list <- read.table("candidate_1.5LOD_8.txt",header=FALSE)

candidate_list <- as.character(candidate_list$V1) 
length(geneID2GO)
length(candidate_list)


#Tell TopGO where these interesting genes appear in the 'geneUniverse' vector, which are all our genes in our peaks of interest here 
geneNames <- names(geneID2GO)
geneList=factor(as.integer(geneNames %in% candidate_list))
names(geneList) <- geneNames
str(geneList)

#Put data in an R object of type 'topGOdata' which will contain the list of genes of interest, the GO annotations, and the GO hierarchy itself. Use MF to denote "Molecular Function" for ontology. 
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata


#List of genes of interest can be accessed using the method sigGenes():
#sg <- sigGenes(myGOdata)
#str(sg)
#numSigGenes(myGOdata)

##ENRICHMENT

#Run Fisher's exact test based on gene counts and take GO hierarchy into account using weight01
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher") 
resultFisher

#List the top ten significant results found
allRes <- GenTable(myGOdata, Fis = resultFisher, orderBy = "resultFisher", ranksOf = "Fis", topNodes = 8)
allRes

#Chromosome 8 
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 8, useInfo ='all')
#GO.ID                                        Term Annotated Significant Expected   Fis
#1  GO:0004333                 fumarate hydratase activity         1           1     0.02 0.017
#2  GO:0004864      protein phosphatase inhibitor activity         1           1     0.02 0.017
#3  GO:0004713            protein tyrosine kinase activity       531          17     9.27 0.021
#4  GO:0016972                      thiol oxidase activity         2           1     0.03 0.035
#5  GO:0004109         coproporphyrinogen oxidase activity         2           1     0.03 0.035
#6  GO:0008378              galactosyltransferase activity        17           2     0.30 0.035
#7  GO:0005319                  lipid transporter activity         3           1     0.05 0.051
#8  GO:0042393                             histone binding         3           1     0.05 0.051

showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 8, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 8, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)


### Id genes annotated with significant Go terms

myterms =allRes$GO.ID
mygenes = genesInTerm(myGOdata, myterms)

var=c()
for (i in 1:length(myterms))
{
  myterm=myterms[i]
  mygenesforterm= mygenes[myterm][[1]]
  mygenesforterm=paste(mygenesforterm, collapse=',')
  var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}


##################


#Chromosome 21 
candidate_list <- read.table("candidate_1.5LOD_21.txt",header=FALSE)

candidate_list <- as.character(candidate_list$V1) 
length(geneID2GO)
length(candidate_list)


#Tell TopGO where these interesting genes appear in the 'geneUniverse' vector, which are all our genes in our peaks of interest here 
geneNames <- names(geneID2GO)
geneList=factor(as.integer(geneNames %in% candidate_list))
names(geneList) <- geneNames
str(geneList)

#Put data in an R object of type 'topGOdata' which will contain the list of genes of interest, the GO annotations, and the GO hierarchy itself. Use MF to denote "Molecular Function" for ontology. 
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata


#List of genes of interest can be accessed using the method sigGenes():
sg <- sigGenes(myGOdata)
#str(sg)
#numSigGenes(myGOdata)

##ENRICHMENT

#Run Fisher's exact test based on gene counts and take GO hierarchy into account using weight01
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher") 
resultFisher

#List the top ten significant results found
allRes <- GenTable(myGOdata, Fis = resultFisher, orderBy = "resultFisher", ranksOf = "Fis", topNodes = 14)
allRes

#Chrome 21
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 14, useInfo ='all')
#GO.ID                                        Term Annotated Significant Expected   Fis
#1  GO:0004180                   carboxypeptidase activity        11           2     0.16 0.011
#2  GO:0004672                     protein kinase activity       547           9     8.11 0.013
#3  GO:0004788           thiamine diphosphokinase activity         1           1     0.01 0.015
#4  GO:0050333 thiamine triphosphate phosphatase activi...         1           1     0.01 0.015
#5  GO:0046966    nuclear thyroid hormone receptor binding         1           1     0.01 0.015
#6  GO:0004799               thymidylate synthase activity         1           1     0.01 0.015
#7  GO:0003834 beta-carotene 15,15'-dioxygenase activit...         1           1     0.01 0.015
#8  GO:0004820                glycine-tRNA ligase activity         1           1     0.01 0.015
#9  GO:0003904 deoxyribodipyrimidine photo-lyase activi...         2           1     0.03 0.029
#10 GO:0004613 phosphoenolpyruvate carboxykinase (GTP) ...         2           1     0.03 0.029
#11 GO:0046983               protein dimerization activity        98           5     1.45 0.033
#12 GO:0003908 methylated-DNA-[protein]-cysteine S-meth...         3           1     0.04 0.044
#13 GO:0004420 hydroxymethylglutaryl-CoA reductase (NAD...         3           1     0.04 0.044
#14 GO:0008484           sulfuric ester hydrolase activity         3           1     0.04 0.044

showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 14, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 14, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)



#Chromosome 21
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 14, useInfo ='all')
#GO.ID                                        Term Annotated Significant Expected   Fis


#all loci 
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 8, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 8, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#1  GO:1990239                     steroid hormone binding         2           2     0.11 0.0032
#2  GO:0030284          nuclear estrogen receptor activity         3           2     0.17 0.0091
#3  GO:0005524                                 ATP binding       123          13     6.92 0.0204
#4  GO:0008081       phosphoric diester hydrolase activity        50           6     2.81 0.0278
#5  GO:0004180                   carboxypeptidase activity        11           3     0.62 0.0281
#6  GO:0004713            protein tyrosine kinase activity       531          41    29.85 0.0346
#7  GO:0003755 peptidyl-prolyl cis-trans isomerase acti...        15           3     0.84 0.0485
#8  GO:0004672                     protein kinase activity       547          44    30.75 0.0521





