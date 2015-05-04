library(igraph)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#====================Load Data====================
patents.df = read.csv("data//patents//apat63_99.txt")
patents.df = patents.df[,c("PATENT", "GYEAR", "COUNTRY", "POSTATE", "NCLASS", "CAT")]
patents.df = patents.df[patents.df$GYEAR>=1995,]
patents.df = patents.df[patents.df$COUNTRY=='US',c("PATENT", "GYEAR", "POSTATE", "NCLASS", "CAT")]

citations.df = read.csv("data//patents//cite75_99.txt")

data.frame(table(patents.df$GYEAR))


nrow(citations.df[(citations.df$CITING>=5590420) & (citations.df$CITED>=5590420),])

citations.df.del = citations.df[(citations.df$CITING %in% patents.df$PATENT) & (citations.df$CITED%in% patents.df$PATENT),]


G <- graph.data.frame(citations.df.del, directed=T)
deg <- degree(G)

tkplot(G)
