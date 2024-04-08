# GO and KEGG enrichment with clusterProfiler
library(clusterProfiler)
library(ggplot2)

# Input file
stat<-read.table("pseudobulk/Mono_stat.txt", header = T, sep = "\t")

genecut<-0.05	# DEG cut off
enrichcut<-0.05	# enrichment cut off

# Significant DEGs (can distinguish down- and up-regulated DEGs)
sig<-subset(stat, padj < genecut)
sig<-subset(stat, padj < genecut & log2FoldChange > 0)

# Gene symbol to ID
symbols<-rownames(sig)

map<-bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id<-map$ENTREZID

rst<-data.frame()

##### PART I. GO enrichment #####
# GO MF enrichment
ego.mf<-enrichGO(gene = id, OrgDb = "org.Hs.eg.db", ont = "MF", pvalueCutoff = enrichcut)
ego.mf<-simplify(ego.mf)								# Remove redundancy
ego.mf<-setReadable(ego.mf, 'org.Hs.eg.db', 'ENTREZID')	# convert ID to symbol
ego.mf<-as.data.frame(ego.mf)

if (nrow(ego.mf) > 0) rst<-rbind(rst, data.frame(Function = "GO_MF", ego.mf))  

# GO BP enrichment
ego.bp<-enrichGO(gene = id, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = enrichcut)
ego.bp<-simplify(ego.bp)
ego.bp<-setReadable(ego.bp, 'org.Hs.eg.db', 'ENTREZID')
ego.bp<-as.data.frame(ego.bp)

if (nrow(ego.bp) > 0) rst<-rbind(rst, data.frame(Function = "GO_BP", ego.bp))  

##### PART II. KEGG enrichment #####
ekegg<-enrichKEGG(gene = id, organism = 'hsa', pvalueCutoff = enrichcut)
ekegg<-setReadable(ekegg, 'org.Hs.eg.db', 'ENTREZID')
ekegg<-as.data.frame(ekegg)
if (nrow(ekegg) > 0) rst<-rbind(rst, data.frame(Function = "KEGG", ekegg))

##### PART III. Hallmark gene sets enrichment ##### 
# Read hallmark gmt file
gmt<-read.gmt("h.all.v2023.1.Hs.symbols.gmt")
egmt<-enricher(symbols, pvalueCutoff = enrichcut, TERM2GENE = gmt)
egmt<-as.data.frame(egmt)
if (nrow(egmt) > 0) rst<-rbind(rst, data.frame(Function = "Hallmark", egmt))

write.table(x = rst, file = "Mono_up_enrich.txt", sep = "\t", quote = F, row.names = F)
