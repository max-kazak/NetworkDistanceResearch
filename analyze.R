library(igraph)
library(ggplot2)
library(reshape2)

source('metrics.R')

#Load Les Miserables graph
load("miserables.Rdata")

#view graph
V(miserables)
plot(miserables)

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
V(miserables)$color = miserables.clust
miserables$layout <- layout.fruchterman.reingold
plot(miserables)
#tkplot(miserables)

#=======================classification=======================================
source('classification.R')

class = classify.multiple(times = 100, graph = miserables, distance = logforest_dist, alpha = 0.01)

#plot
V(miserables)$color = class
plot(miserables)

#=====================modularity and error alpha-resistance==========================
alphas <- seq(0.01, 0.91, by=0.01)
dist.vect <- c(plainwalk_dist, walk_dist , plainforest_dist, logforest_dist, communicability_dist, logcommunicability_dist)
names(dist.vect) <- c("plain_walk", "walk", "plain_forest", "log_forest", "communicability", "log_communicability")

mods.df = data.frame(alpha = alphas)
acc.df = data.frame(alpha = alphas)

#calculate modularity for each alpha in alphas
for(i in 1:length(dist.vect)) {
    mods = vector()
    acc = vector()
    for(a in alphas) {
        class = classify.multiple(times = 100, graph = miserables, distance = dist.vect[[i]], alpha = a)
        mods = c(mods, modularity(miserables, class))
        acc = c(acc, sum(class==miserables.clust)/length(class))
    }
    mods.df[names(dist.vect[i])] = mods
    acc.df[names(dist.vect[i])] = acc
}
mods.df
acc.df

#view dependency of alpha and modularity
mods.melted.df <- melt(mods.df, id=c("alpha"))
g.mods<-ggplot(mods.melted.df) + 
    geom_point(aes(alpha, value, colour=variable)) +
    geom_smooth(aes(alpha, value, colour=variable)) +
    geom_line(aes(alpha, value, colour=variable))
plot(g.mods)

#view dependency of alpha and error
acc.melted.df <- melt(acc.df, id=c("alpha"))
g.acc<-ggplot(acc.melted.df) + 
                 geom_point(aes(alpha, value, colour=variable)) +
                 geom_smooth(aes(alpha, value, colour=variable)) +
                 geom_line(aes(alpha, value, colour=variable))
plot(g.acc)


