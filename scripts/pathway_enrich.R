library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

print("Start loading...")
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
# filename <- "./results/SE_MATS_JC.csv"

event_name <- unlist(strsplit(
  unlist(strsplit(filename, split="/"))[3], split="_"))[1]

mat <- read.csv(filename, row.names = "gene_id")
print("data loaded!")

geneUniverse <- unlist(rownames(mat))

mask <- mat$p_value < 0.05
deGenes <- unlist(rownames(mat[mask, ]))

#=========================== Enrichment with GO ===========================

# readable must be False if using "SYMBOL"
print("Start enrichment with GO...")
ans.go.bp <- enrichGO(gene = deGenes, ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   keyType = "SYMBOL",
                   universe = geneUniverse,
                   readable=FALSE,
                   pvalueCutoff = 0.05)
ans.go.cc <- enrichGO(gene = deGenes, ont = "CC",
                      OrgDb ="org.Hs.eg.db",
                      keyType = "SYMBOL",
                      universe = geneUniverse,
                      readable=FALSE,
                      pvalueCutoff = 0.05)
ans.go.mf <- enrichGO(gene = deGenes, ont = "MF",
                      OrgDb ="org.Hs.eg.db",
                      keyType = "SYMBOL",
                      universe = geneUniverse,
                      readable=FALSE,
                      pvalueCutoff = 0.05)
print("Enrichment with Go Done!")

tab.go.bp <- as.data.frame(ans.go.bp)
tab.go.bp <- subset(tab.go.bp, Count>5)
#p1 <- barplot(ans.go.bp, showCategory=10, title = "AD vs CTRL (GO-BP)")

tab.go.cc <- as.data.frame(ans.go.cc)
tab.go.cc <- subset(tab.go.cc, Count>5)
#p2 <- barplot(ans.go.cc, showCategory=10, title = "AD vs CTRL (GO-CC)")

tab.go.mf <- as.data.frame(ans.go.mf)
tab.go.mf <- subset(tab.go.mf, Count>5)
#p3 <- barplot(ans.go.mf, showCategory=10, title = "AD vs CTRL (GO-MF)")


if(nrow(tab.go.bp) > 30){
  ans.go.bp <- simplify(ans.go.bp, 
                               cutoff=0.5,
                               by="p.adjust",
                               select_fun=min)
  tab.ans.go.bp <- as.data.frame(ans.go.bp)
}

enrichment_result_dir = paste(event_name, "event_enrichment", sep="_")
dir.create(enrichment_result_dir)
if (nrow(tab.go.bp) > 0){
  write.csv(tab.go.bp,
            file.path(enrichment_result_dir, paste(event_name, "BP.csv", sep="_")),
            row.names = FALSE)
  temp <- pairwise_termsim(ans.go.bp)
  pdf(file.path(enrichment_result_dir, paste("BP", event_name, "EMA.pdf", sep="_")))
  plot(emapplot(temp))
  dev.off()
  pdf(file.path(enrichment_result_dir, paste("BP", event_name, "Dot.pdf", sep="_")))
  plot(dotplot(ans.go.bp, showCategory=15) + ggtitle("GO-BP"))
  dev.off()
}

if (nrow(tab.go.cc) > 0){
  write.csv(tab.go.cc,
            file.path(enrichment_result_dir, paste(event_name, "CC.csv", sep="_")),
            row.names = FALSE)
  temp <- pairwise_termsim(ans.go.cc)
  pdf(file.path(enrichment_result_dir, paste("CC", event_name, "EMA.pdf", sep="_")))
  plot(emapplot(temp))
  dev.off()
  pdf(file.path(enrichment_result_dir, paste("CC", event_name, "Dot.pdf", sep="_")))
  plot(dotplot(ans.go.cc, showCategory=15) + ggtitle("GO-CC"))
  dev.off()
}

if (nrow(tab.go.mf) > 0){
  write.csv(tab.go.mf,
            file.path(enrichment_result_dir, paste(event_name, "MF.csv", sep="_")),
            row.names = FALSE)
  temp <- pairwise_termsim(ans.go.mf)
  pdf(file.path(enrichment_result_dir, paste("MF", event_name, "EMA.pdf", sep="_")))
  plot(emapplot(temp))
  dev.off()
  pdf(file.path(enrichment_result_dir, paste("MF", event_name, "Dot.pdf", sep="_")))
  plot(dotplot(ans.go.mf, showCategory=15) + ggtitle("GO-MF"))
  dev.off()
}

# #============================================================================
# #=========================== Enrichment with KEGG ===========================
# print("Start enrichment with KEGG...")
# entrezGenes <- unlist(mget(deGenes, envir=org.Hs.egSYMBOL2EG,
#                        ifnotfound = NA))

# ans.kegg <- enrichKEGG(gene = entrezGenes,
#                        organism = "hsa",
#                        universe = geneUniverse,
#                        pvalueCutoff = 0.05,
#                        minGSSize = 0, maxGSSize = 500000,
#                        use_internal_data = TRUE)

# tab.kegg <- as.data.frame(ans.kegg)
# tab.kegg<- subset(tab.kegg, Count>5)

#p4 <- barplot(ans.kegg, showCategory=10, title = "AD vs CTRL (KEGG)")
