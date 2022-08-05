
#set working directory
setwd("C:/Users/Susan Coort/Documents/Papers/VitaminD_Chapter/Pathway_analysis")

#load required libraries
library(dplyr)
library(RColorBrewer)
library(org.Hs.eg.db)
library(tidyverse)
library(EnhancedVolcano)
library(VennDiagram)
library(rWikiPathways)
library(data.table)
library(RCy3)

# Differential expression data visualization
# We will use a publicly available dataset, which identified two different 
# subtypes of IBD - CD and UC. DEG analysis was performed before then we will use the results file from prevous workflow
# First, let's import the data and use a volcano plot to visualize the result 
# of the differential gene expression analysis result, and use a Venn diagram 
# to study how many differentially expressed genes are shared between the subtypes.

#DC read, the datasets may change depends on our observation (ileum or rectum)
dataset.DC <- read.delim("DC_comp_Group_immature_vitD-Group_immature_control.txt")
#mono read
dataset.mono <- read.delim("monocyte_comp_Group__vitaminD-Group_control.txt")
#THP1 read
dataset.thp1 <- read.delim("THP1_comp_Group_Uninfected_treated_VitD-Group_Uninfected_untreated.txt")

#get SYMBOL of each gene symbols for each disease type
hs <- org.Hs.eg.db
Symbol.DC <- AnnotationDbi::select(hs, 
            keys = dataset.DC$ENSG_ID,
            columns = c("SYMBOL", "ENTREZID"),
            keytype = "ENSEMBL")


Symbol.mono <- AnnotationDbi::select(hs, 
                                   keys = dataset.mono$ENSG_ID,
                                   columns = c("SYMBOL", "ENTREZID"),
                                   keytype = "ENSEMBL")

Symbol.thp1 <- AnnotationDbi::select(hs, 
                                     keys = dataset.thp1$ENSG_ID,
                                     columns = c("SYMBOL", "ENTREZID"),
                                     keytype = "ENSEMBL")


#filter out double gene symbols
Symbol.DC<- Symbol.DC%>% distinct(Symbol.DC$ENSEMBL, .keep_all = TRUE)
Symbol.mono<- Symbol.mono%>% distinct(Symbol.mono$ENSEMBL, .keep_all = TRUE)
Symbol.thp1<- Symbol.thp1%>% distinct(Symbol.thp1$ENSEMBL, .keep_all = TRUE)

#get back up data before starting
#dataset.b  <- dataset.DC
#dataset.b2 <- dataset.mono

# add SYmbols to each dataset
dataset.DC <- cbind(Symbol.DC$SYMBOL, Symbol.DC$ENTREZID,dataset.DC)
dataset.mono <- cbind(Symbol.mono$SYMBOL, Symbol.mono$ENTREZID,dataset.mono)
dataset.thp1 <- cbind(Symbol.thp1$SYMBOL, Symbol.thp1$ENTREZID, dataset.thp1)

#change column names
colnames(dataset.DC)[1] <- "SYMBOL"
colnames(dataset.DC)[2] <- "ENTREZID"
colnames(dataset.mono)[1] <- "SYMBOL"
colnames(dataset.mono)[2] <- "ENTREZID"
colnames(dataset.thp1)[1] <- "SYMBOL"
colnames(dataset.thp1)[2] <- "ENTREZID"



# filter genes without Entrez Gene identifier
#dataset.CD <- dataset.CD %>% tidyr::drop_na(ENTREZ.ID)
#dataset.UC <- dataset.UC %>% tidyr::drop_na(ENTREZ.ID)
#remove some unused columns
#dataset.CD <- subset( dataset.CD, select = -c(3,5,6,7,9 ) )
#dataset.UC <- subset( dataset.UC, select = -c(3,5,6,7,9 ) )

#output folder should be created beforehand run the code
png('output/volcanoplot_DC.png')
EnhancedVolcano(dataset.DC , title = "DC", lab = dataset.DC$SYMBOL, 
                labSize = 3, x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

png('output/volcanoplot_mono.png')
EnhancedVolcano(dataset.mono, title = "Monocytes", lab = dataset.mono$SYMBOL, 
                labSize = 3, x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

png('output/volcanoplot_THP1.png')
EnhancedVolcano(dataset.thp1, title = "THP1", lab = dataset.thp1$SYMBOL, 
                labSize = 3, x = 'logFC', y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

#list of all deg from two disease types 
deg.DC <-unique(dataset.DC[!is.na(dataset.DC$adj.P.Val) & dataset.DC$adj.P.Val < 0.05 & abs(dataset.DC$logFC) > 0.26 & dataset.DC$AveExpr > 3.86, c(1,2,3,4)])
DC.up   <-unique(dataset.DC[!is.na(dataset.DC$adj.P.Val) & dataset.DC$adj.P.Val < 0.05 & dataset.DC$logFC > 0.26 & dataset.DC$AveExpr > 3.86, c(1,2,3,4)])
DC.down <-unique(dataset.DC[!is.na(dataset.DC$adj.P.Val) & dataset.DC$adj.P.Val < 0.05 & dataset.DC$logFC < -0.26 & dataset.DC$AveExpr > 3.86, c(1,2,3,4)])

deg.mono <-unique(dataset.mono[!is.na(dataset.mono$adj.P.Val) & dataset.mono$adj.P.Val < 0.05 & abs(dataset.mono$logFC) > 0.26 & dataset.mono$AveExpr > 4.7, c(1,2,3,4)])
mono.up   <-unique(dataset.mono[!is.na(dataset.mono$adj.P.Val) & dataset.mono$adj.P.Val < 0.05 & dataset.mono$logFC > 0.26 & dataset.mono$AveExpr > 4.7, c(1,2,3,4)])
mono.down <-unique(dataset.mono[!is.na(dataset.mono$adj.P.Val) & dataset.mono$adj.P.Val < 0.05 & dataset.mono$logFC < -0.26 & dataset.mono$AveExpr > 4.7, c(1,2,3,4)])


deg.thp1 <-unique(dataset.thp1[!is.na(dataset.thp1$adj.P.Val) & dataset.thp1$adj.P.Val < 0.05 & abs(dataset.thp1$logFC) > 0.26 & dataset.thp1$AveExpr > 4.25, c(1,2,3,4)])
thp1.up   <-unique(dataset.thp1[!is.na(dataset.thp1$adj.P.Val) & dataset.thp1$adj.P.Val < 0.05 & dataset.thp1$logFC > 0.26 & dataset.thp1$AveExpr > 4.25, c(1,2,3,4)])
thp1.down <-unique(dataset.thp1[!is.na(dataset.thp1$adj.P.Val) & dataset.thp1$adj.P.Val < 0.05 & dataset.thp1$logFC < -0.26 & dataset.thp1$AveExpr > 4.25, c(1,2,3,4)])

venn.diagram(x = list(deg.DC$ENSG_ID, deg.mono$ENSG_ID, deg.thp1$ENSG_ID),
             category.names = c("DC" ,"Mono", "THP-1"),
             filename = 'output/venn_DEGs.png',
             output=FALSE,
             col=c("#440154ff","#440154ff", '#21908dff'),
             fill = c(alpha("#440154ff",0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
             cex = 1.5,
)

## Pathway enrichment analysis
#We will perform pathway enrichment with the gene sets of all pathway models in WikiPathways (human only).
gmt <- "wikipathways-20220310-gmt-Homo_sapiens.gmt"
wp2gene   <- readPathwayGMT(gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
bkgd.genes <- unique(dataset.DC[,c(1,2)])# genes are commonn for both disease type CD and UC

# The clusterProfiler R-package is used to perform overrepresentation analysis (ORA). 
# The function can be easily replaced to use other enrichment methods (GSEA / rSEA / etc). 
# We will run the analysis separately for CD and UC subtype.

##################################FOR DC ##############################################
ewp.DC <- clusterProfiler::enricher(
  deg.DC$ENTREZID,
  universe = bkgd.genes$ENTREZID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.DC.res <- as.data.frame(ewp.DC) 

# number of genes measured in pathways
length(ewp.DC@universe)
# number of DEG in pathways
length(deg.DC$ENTREZID[deg.DC$ENTREZID %in% unique(wp2gene$gene)])
num.pathways.DC <- dim(ewp.DC.res)[1]


# export enrichment result
png('output/DC_barplot.png', width = 1500, height=1000)
ggplot(ewp.DC[1:num.pathways.DC], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#ADD8E6") +
  coord_flip() +
  labs(x="", y="DC DEG gene count", fill="") +
  theme(axis.text=element_text(size=20)) + 
  theme(legend.position="none")
dev.off()
write.table(ewp.DC.res, file="output/DC_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)

# > Interpretation
# - **Q5**: How many pathways are altered in the CD subtype and how do they link to 
#IBD  (expected or unexpected)?

##################################FOR Monocytes ##############################################
ewp.mono <- clusterProfiler::enricher(
  deg.mono$ENTREZID,
  universe = bkgd.genes$ENTREZID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.mono.res <- as.data.frame(ewp.mono) 

# number of genes measured in pathways
length(ewp.mono@universe)
# number of DEG in pathways
length(deg.mono$ENTREZID[deg.mono$ENTREZID %in% unique(wp2gene$gene)])
num.pathways.mono <- dim(ewp.mono.res)[1]

# export enrichment result
png('output/mono_barplot.png', width = 1500, height=1000)
ggplot(ewp.mono[1:num.pathways.mono], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#21908dff") +
  coord_flip() +
  labs(x="", y="Mono DEG gene count", fill="") +
  theme(axis.text=element_text(size=20)) + 
  theme(legend.position="none")
dev.off()
write.table(ewp.mono.res, file="output/Mono_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)



##################################FOR THP1 ##############################################
bkgd.genes <- unique(dataset.thp1[,c(1,2)])


ewp.thp1 <- clusterProfiler::enricher(
  deg.thp1$ENTREZID,
  universe = bkgd.genes$ENTREZID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.thp1.res <- as.data.frame(ewp.thp1) 

# number of genes measured in pathways
length(ewp.thp1@universe)
# number of DEG in pathways
length(deg.thp1$ENTREZID[deg.thp1$ENTREZID %in% unique(wp2gene$gene)])
num.pathways.thp1 <- dim(ewp.thp1.res)[1]

# export enrichment result
png('output/THP1_barplot.png', width = 1500, height=1000)
ggplot(ewp.thp1[1:num.pathways.thp1], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#440154ff") +
  coord_flip() +
  labs(x="", y="THP-1 DEG gene count", fill="") +
  theme(axis.text=element_text(size=15)) + 
  theme(legend.position="none")
dev.off()
write.table(ewp.thp1.res, file="output/THP1_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)




#############Pathway visualization##############
# The pathways can then be visualized with the gene expression data as shown with the 
# "Overview of proinflammatory and profibrotic mediators" (WP5095) pathway from WikiPathways. 
# The pathway was altered in both subtypes. 
# interest and visualize the data on that pathway. 

#before starting we should merge data table of CD and UC disease 
data.immuno <- merge(dataset.DC, dataset.mono, by = "SYMBOL")

colnames(data.immuno)[2] <- "ENTREZID"
colnames(data.immuno)[4] <- "logFC_DC"
colnames(data.immuno)[8] <- "pvalue_DC"
colnames(data.immuno)[9] <- "adj,pvalue_DC"
colnames(data.immuno)[13] <- "logFC_Mono"
colnames(data.immuno)[17] <- "pvalue_Mono"
colnames(data.immuno)[18] <- "adj.pvalue_Mono"


RCy3::cytoscapePing()
RCy3::installApp(c("wikipathways","CyTargetLinker"))

#glycolysis and glyconeogensis pathway
RCy3::commandsRun('wikipathways import-as-pathway id=WP534') 
toggleGraphicsDetails()
RCy3::mapTableColumn("Ensembl", "Homo sapiens", "Ensembl", "Entrez Gene")
loadTableData(data.immuno, data.key.column = "ENSG_ID.x", table.key.column = "Ensembl")

RCy3::installApp("enhancedGraphics")
RCy3::copyVisualStyle("WikiPathways", "my_style_heatmap")
RCy3::setNodeCustomHeatMapChart(c("logFC_DC","logFC_Mono"), slot = 2, style.name = "my_style_heatmap", 
                                colors = c("#CC3300","#FFFFFF","#6699FF","00000000"))

RCy3::setVisualStyle("my_style_heatmap")

######################################## Oxidative phosphorylation ##########################################
RCy3::commandsRun('wikipathways import-as-pathway id=WP111') 
toggleGraphicsDetails()
loadTableData(data.immuno, data.key.column = "ENSG_ID.x", table.key.column = "Ensembl")

RCy3::installApp("enhancedGraphics")
RCy3::copyVisualStyle("WikiPathways", "my_style_heatmap")
RCy3::setNodeCustomHeatMapChart(c("logFC_DC","logFC_Mono"), slot = 2, style.name = "my_style_heatmap", 
                                colors = c("#CC3300","#FFFFFF","#6699FF","00000000"))

RCy3::setVisualStyle("my_style_heatmap")



#Saving output
png.file <- file.path("output/PathwayVisualization.png")
exportImage(png.file,'PNG', zoom = 500)
cys.file <- file.path(getwd(),"output/PathwayVisualization.cys")
saveSession(cys.file)
#comment following line if you want to manipulate the visualization in Cytoscape
RCy3::closeSession(save.before.closing = F)


write.table(data.immuno, file="output/merged_data.txt", sep="\t", quote=FALSE, row.names = FALSE)




##########################Pathway overlap visualization#################################
# There is often crosstalk and overlap between pathways enriched in gene expression analyses. 
# The following step visualizes the overlap between the enriched pathways in a pathway-gene network. 
# The genes not present in any pathway are included in the visualization but can be filtered 
#in a follow-up step if preferred. 

pwy <- unique(ewp.DC.res[,c(1,2)])#enriched pathways in CD 
colnames(pwy) <- c("id","label")
pwy$type <- 'pathway'
edges <- wpid2gene[wpid2gene$wpid %in% pwy$id,]
colnames(edges) <- c("source", "target")
genes <- unique(deg.DC)
colnames(genes) <- c("id","label")
genes <- transform(genes, id = as.character(id))
genes$type <- 'gene'
edges <- unique(edges[edges$target %in% genes$id,])
genes <- genes[genes$id %in% edges$target,]
nodes <- dplyr::bind_rows(genes, pwy)
rownames(nodes) <- NULL
createNetworkFromDataFrames(nodes=nodes,edges=edges,title="Pathway-Gene-Associations", collection="PathwayGeneCrosstalk")
loadTableData(data.immuno, data.key.column = "ENTREZID", table.key.column = "id")

# Visual style
RCy3::copyVisualStyle("default","wp.vis")
RCy3::setNodeLabelMapping("label", style.name="wp.vis")
RCy3::lockNodeDimensions(TRUE, style.name="wp.vis")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="wp.vis")
RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "wp.vis")
data.values<-c(-1,0,1) 
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("logFC_CD", data.values, node.colors, default.color = "#99FF99", style.name = "wp.vis")
RCy3::setVisualStyle("wp.vis")
RCy3::toggleGraphicsDetails()

# Saving output
svg.file <- file.path(getwd(), "output/PathwayCrosstalk.svg")
exportImage(svg.file,'SVG')
png.file <- file.path(getwd(), "output/PathwayCrosstalk.png")
exportImage(png.file,'PNG', zoom = 500)
cys.file <- file.path(getwd(), "output/PathwayCrosstalk.cys")
saveSession(cys.file) 
#comment following line if you want to manipulate the visualization in Cytoscape
#RCy3::closeSession(save.before.closing = F)

##############################Drug target information#####################################
# Next, we will add information about known drug-target interactions for the genes 
# in the affected pathways using information from DrugBank using the CyTargetLinker app.
# We will show this for the UC subtype. 

RCy3::cytoscapePing()
installApp('CyTargetLinker') 

pwy <- unique(ewp.CD.res[,c(1,2)])
colnames(pwy) <- c("id","label")
pwy$type <- 'pathway'
edges <- wpid2gene[wpid2gene$wpid %in% pwy$id,]
colnames(edges) <- c("source", "target")
genes <- unique(deg.CD)
colnames(genes) <- c("id","label")
genes <- transform(genes, id = as.character(id))
genes$type <- 'gene'
edges <- unique(edges[edges$target %in% genes$id,])
genes <- genes[genes$id %in% edges$target,]
nodes <- dplyr::bind_rows(genes, pwy)
rownames(nodes) <- NULL
createNetworkFromDataFrames(nodes=nodes,edges=edges,title="Pathway-Gene-Associations", collection="PathwayGeneCrosstalk")
loadTableData(data.immuno, data.key.column = "ENTREZID", table.key.column = "id")

# Visual style
RCy3::copyVisualStyle("default","wp.vis")
RCy3::setNodeLabelMapping("label", style.name="wp.vis")
RCy3::lockNodeDimensions(TRUE, style.name="wp.vis")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="wp.vis")
RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "wp.vis")
data.values<-c(-1,0,1) 
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("logFC_UC", data.values, node.colors, default.color = "#99FF99", style.name = "wp.vis")
RCy3::setVisualStyle("wp.vis")
RCy3::toggleGraphicsDetails()

drugbank <- file.path(getwd(), "data/drugbank-5.1.0.xgmml")

# run CyTargetLinker
commandsRun(paste0('cytargetlinker extend idAttribute="id" linkSetFiles="', drugbank, '"'))
commandsRun('cytargetlinker applyLayout network="current"')
RCy3::setVisualStyle("wp.vis")

#let's change the visualization of the drugs in the network using the ByPass option
selected <- RCy3::selectNodes(nodes="drug", by.col = "CTL.Type")
RCy3::setNodeShapeBypass(node.names = selected$nodes, new.shapes = "Triangle")
RCy3::setNodeColorBypass(node.names = selected$nodes, "#FFFFCE")
RCy3::setNodeBorderColorBypass(node.names = selected$nodes, "#000000")
RCy3::setNodeBorderWidthBypass(node.names = selected$nodes, 4)
RCy3::clearSelection()
RCy3::toggleGraphicsDetails()

png.file <- file.path(getwd(), "output/drug_target.png")
exportImage(png.file,'PNG', zoom = 500)






