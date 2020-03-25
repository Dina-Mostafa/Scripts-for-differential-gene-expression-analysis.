
#prepare two files with DE genes (TR) and DE proteins (MS) with ensemble gene ID and logFC only
#Upload files 
MS<- read.csv("MS.csv", stringsAsFactors = FALSE)
TR<- read.csv("TR.csv", stringsAsFactors = FALSE)

MS2 <- subset(MS, MS$Gene_ID %in% TR$Gene_ID)
TR2 <- subset(TR, TR$Gene_ID %in% MS$Gene_ID)

all_MS_TR_common<-merge(MS2,TR2, by="Gene_ID")
write.csv(all_MS_TR_common, "All_MS_TR_common.csv", row.names = FALSE, quote= FALSE)

#Correlation analysis
install.packages("ggpubr")
library(ggpubr)

ggscatter(all_MS_TR_common, x = "MS_logFC", y = "TR_logFC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Proteomic log2FC", ylab = "Transcriptomic log2FC")








