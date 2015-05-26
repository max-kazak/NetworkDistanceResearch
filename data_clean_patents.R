library(igraph)
library(ggplot2)

#====================Load Data====================
df.patents = read.csv("data//patents//apat63_99.txt")
df.patents = df.patents[,c("PATENT", "GYEAR", "COUNTRY", "POSTATE", "NCLASS", "CAT")]
df.citations = read.csv("data//patents//cite75_99.txt")

#===========================filter by year===================
df.patents.del = df.patents[df.patents$COUNTRY=='GB',c("PATENT", "GYEAR", "NCLASS", "CAT")]

data.frame(table(df.patents.del$GYEAR))
df.patents.del = df.patents.del[df.patents.del$GYEAR<=1995,]
df.citations.del = df.citations[(df.citations$CITING %in% df.patents.del$PATENT) & (df.citations$CITED%in% df.patents.del$PATENT),]

G <- graph.data.frame(df.citations.del, directed=F)

#===========================clusters=========================
patents.clust = clusters(G)

head(sort(patents.clust$csize, decreasing = T))
v.patents.clust = as.numeric(V(G)$name)[patents.clust$membership==which.max(patents.clust$csize)]

df.patents.maxclust = df.patents.del[df.patents.del$PATENT %in% v.patents.clust,]
df.citations.maxclust = df.citations.del[(df.citations.del$CITING %in% df.patents.maxclust$PATENT) & 
                                             (df.citations.del$CITED %in% df.patents.maxclust$PATENT),]
G.maxclust <- graph.data.frame(df.citations.maxclust, directed=F)

#========================degree=============
G.deg <- delete.vertices(G.maxclust, V(G.maxclust)[degree(G.maxclust)<=1])

#===========================final======================
v.patents.deg <- as.numeric(V(G.deg)$name)
df.patents.fin = df.patents.maxclust[df.patents.maxclust$PATENT %in% v.patents.deg,]
df.citations.fin = df.citations.maxclust[(df.citations.maxclust$CITING %in% df.patents.fin$PATENT) & 
                                             (df.citations.maxclust$CITED %in% df.patents.fin$PATENT),]
G.fin <- graph.data.frame(df.citations.fin, directed=F)

plot(G.fin, vertex.label=NA, vertex.size=5)

cat = sapply(as.numeric(V(G.fin)$name), function(pat) df.patents.fin[df.patents.fin$PATENT==pat,]$CAT)

V(G.fin)$color = cat;


