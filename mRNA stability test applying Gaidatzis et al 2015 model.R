# load edgeR library
library(edgeR)

# input files and parameters
exsFile <- "exon_reads.txt"  #exons raw counts file

conditions <- c("Control","Control","Control","KO", "KO", "KO") # correspond to the columns in exsFile and insFile; the first condition will be the reference

# read in count tables and remove the first column (ensembl ID)
cntEx <- read.delim(exsFile)[,-1]

# combine exon and intron raw counts for edgeR containing genes with sufficient counts in both exonic and intronic levels
raw_counts <- read.csv("Cnt.csv", stringsAsFactors=FALSE)  # Cnt is combined exons and introns file 240319 keeping only normalized exons and introns and ensembl IDs
Cnt <- raw_counts[ , -1]
rownames(Cnt)<-raw_counts[ , 1] #gene Ensembl IDs
View(Cnt)

# edgeR workflow
factorRegion    <- factor(rep(c("ex","in"),each=ncol(cntEx)), levels=c("in", "ex")) # define experimental factor exon/intron 
factorCondition <- factor(rep(conditions,2), levels=unique(conditions)) # define experimental factor 'conditions'

y <- DGEList(counts=Cnt, genes=data.frame(EnsemblID=rownames(Cnt))) # define DGEList object
y <- calcNormFactors(y) # determine normalization factors
design <- model.matrix(~ factorRegion * factorCondition) # design matrix with interaction term
rownames(design) <- colnames(Cnt)

y <- estimateDisp(y, design) # estimate dispersion
fit <- glmFit(y, design) # fit generalized linear model
lrt <- glmLRT(fit) # calculate likelihood-ratio between full and reduced models
#final table with significance level for each gene 
tt <- topTags(lrt, n=nrow(y))
head(tt$table)
write.csv(tt$table, file = "tt.csv", row.names = FALSE, quote = FALSE )

stablitytable <- read.csv("tt.csv", stringsAsFactors=FALSE)
#to add gene symbols to table
gene_ref <- read.csv("gene_ref.csv", stringsAsFactors=FALSE)

final<- merge(stablitytable, gene_ref, by = "Gene_ID")
write.csv(final, file = "Final stability.csv", row.names = FALSE, quote = FALSE )
