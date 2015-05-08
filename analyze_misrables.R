library(igraph)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#prepare for parallel computation
library(doParallel)
library(foreach)
numWorkers=4

source('metrics.R')

#Load Les Miserables graph
load("data/miserables/miserables.Rdata")
names <- read.table("data/miserables/names.txt", stringsAsFactors=F)$V1
V(miserables)$name <- names

#add fixed layout
mylayout <- layout.fruchterman.reingold(miserables)
miserables$layout <- mylayout

#choose palette
palette(brewer.pal(7, "Dark2"))
#display.brewer.all()

#view graph
V(miserables)
plot(miserables, vertex.shape="none")

# let's see if we have communities here using the 
# Grivan-Newman algorithm
# 1st we calculate the edge betweenness, merges, etc...
ebc <- edge.betweenness.community(miserables, directed=F)

# Now we have the merges/splits and we need to calculate the modularity
# for each merge for this we'll use a function that for each edge
# removed will create a second graph, check for its membership and use
# that membership to calculate the modularity
mods <- sapply(0:ecount(miserables), function(i){
    g <- delete.edges(miserables, ebc$removed.edges[seq(length=i)])
    cl <- clusters(g)$membership
    modularity(miserables,cl)
})
# we can now plot all modularities
plot(mods, pch=20)

#save clusterization with max modularity
miserables.separated = delete.edges(miserables, ebc$removed.edges[seq(length=which.max(mods)-1)])
miserables.clust = clusters(miserables.separated)$membership
#clusters in df
df.clust = data.frame(V(miserables)$name,miserables.clust)
colnames(df.clust) = c("V","clust")

#view graph colored by cluster
V(miserables)$label.color = miserables.clust
plot(miserables, vertex.shape="none")
#tkplot(miserables)

#=======================classification=======================================
source('classification.R')

l.marked.mult = select.marked.multiple(df.clust)

class = classify.multiple(graph = miserables, marked.mult = l.marked.mult, distance = logforest_dist, alpha = 0.01)

#plot
V(miserables)$label.color = class
plot(miserables, vertex.shape="none")

#=====================modularity and error alpha-resistance==========================
alphas <- seq(0.001, 1.0, by=0.001)
dist.vect <- c(plainwalk_dist, walk_dist , plainforest_dist, logforest_dist, communicability_dist, logcommunicability_dist)
names(dist.vect) <- c("plain_walk", "walk", "plain_forest", "log_forest", "communicability", "log_communicability")

cl = makeCluster(numWorkers)
registerDoParallel(cl)
strt <- Sys.time()
#calculate modularity for each alpha in alphas
lres = foreach(i = 1:length(dist.vect), .packages="igraph") %dopar% {
        source("metrics.R")
        mods = vector()
        acc = vector()
        for(a in alphas) {
            class = classify.multiple(graph = miserables, marked.mult = l.marked.mult, distance = dist.vect[[i]], alpha = a)
            mods = c(mods, modularity(miserables, class))
            acc = c(acc, sum(class==miserables.clust)/length(class))
        }
        data.frame(mods, acc)
}
print(Sys.time()-strt)
stopCluster(cl)

mods.df = data.frame(alpha = alphas)
acc.df = data.frame(alpha = alphas)
for(i in 1:length(dist.vect)) {
    mods.df[names(dist.vect[i])] = lres[[i]]$mods
    acc.df[names(dist.vect[i])] = lres[[i]]$acc
}

#==============print results============
folder = as.character(Sys.Date())
dir.create(folder, showWarnings = FALSE)

#view dependency of alpha and modularity
mods.melted.df <- melt(mods.df, id=c("alpha"))
g.mods<-ggplot(mods.melted.df) + 
    geom_point(aes(alpha, value, colour=variable)) +
    #geom_smooth(aes(alpha, value, colour=variable)) +
    geom_line(aes(alpha, value, colour=variable)) +
    ggtitle("Модульность")
png(filename=file.path(folder,'miserables_mod.png'))
plot(g.mods)
dev.off()

#view dependency of alpha and error
acc.melted.df <- melt(acc.df, id=c("alpha"))
g.acc<-ggplot(acc.melted.df) + 
                 geom_point(aes(alpha, value, colour=variable)) +
                 #geom_smooth(aes(alpha, value, colour=variable)) +
                 geom_line(aes(alpha, value, colour=variable))+
                 ggtitle("Процент Ошибки")
png(filename=file.path(folder,'miserables_error.png'))
plot(g.acc)
dev.off()



