library(dplyr)
library(tidyr)
GWAS.clumped <- read.table("clumped.formatted", header=T, sep="\t")
Genetic_map <- read.table("clumpedSNP.map", header=F, sep="\t"); colnames(Genetic_map) <- c("CHR","SNP","cM","POS")
GWAS.clumped.cM <- left_join(GWAS.clumped, Genetic_map, by = c("CHR" = "CHR","POS" = "POS"))
GWAS.clumped.cM.NA <- subset(GWAS.clumped.cM, is.na(cM))
cat("leads without cM:", nrow(GWAS.clumped.cM.NA), "\n")
GWAS.clumped.cM.ordered <- GWAS.clumped.cM[with(GWAS.clumped.cM, order(CHR, POS)), ]
GWAS.clumped.cM.ordered$Merged <- NA
IndepedentLoci = 1
GWAS.clumped.cM.ordered[1,"Merged"] = IndepedentLoci
for (i in 2:nrow(GWAS.clumped.cM.ordered)){
  if(GWAS.clumped.cM.ordered[i,"CHR"] == GWAS.clumped.cM.ordered[i-1,"CHR"] & abs(GWAS.clumped.cM.ordered[i,"cM"] - GWAS.clumped.cM.ordered[i-1,"cM"]) < 0.1){
  	GWAS.clumped.cM.ordered[i,"Merged"] = IndepedentLoci
  }else {
    IndepedentLoci = IndepedentLoci + 1
  	GWAS.clumped.cM.ordered[i,"Merged"] = IndepedentLoci
  }
}
GWAS.clumped.cM.ordered.Top1 <- GWAS.clumped.cM.ordered %>% group_by(Merged) %>% top_n(n = -1, wt = P) %>% as.data.frame
write.table(GWAS.clumped.cM.ordered[,c("CHR","POS","SNP.x","P","cM","Merged")], "clumped_leads_merged.tsv", col.names=c("CHR","POS","SNP","P","cM","Merged"), row.names=F, quote=F, sep="\t")
write.table(GWAS.clumped.cM.ordered.Top1, "independent_loci.tsv", row.names=F, quote=F, sep="\t")
cat("n_clumped_leads:", nrow(GWAS.clumped.cM.ordered), "\n")
cat("n_independent_loci:", length(unique(GWAS.clumped.cM.ordered$Merged)), "\n")
cat(sprintf('{"n_clumped_leads": %d, "n_independent_loci": %d}\n', nrow(GWAS.clumped.cM.ordered), length(unique(GWAS.clumped.cM.ordered$Merged))), file="summary.json")
