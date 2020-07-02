library("NeatMap")

#read data file

raw_counts <- read.csv("disallowed_raw_counts.csv", stringsAsFactors=FALSE)
data <- raw_counts[ , -1]
rownames(data)<-raw_counts[ , 1] #gene names
View(data)


#data change to log2 scale

data1 <- log2(data)



#generate PDF and heatmap

#pdf("HM.pdf")

heatmap1(as.matrix(data1), row.order = NULL, column.order = NULL, row.cluster = NULL, column.cluster = NULL, row.labels = NULL, column.label.size = 1, row.label.size = 3,row.normalize=F)+scale_fill_gradient2(low="blue", high="red", midpoint = 0)

dev.off()

browseVignettes("heatmaps")
