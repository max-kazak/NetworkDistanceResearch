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
load("miserables.Rdata")
names <- read.table("names.txt", stringsAsFactors=F)$V1
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
# walktrap algorithm
wc <- walktrap.community(miserables)

#save clusterization with max modularity
miserables.clust = wc$membership
#clusters in df
df.clust = data.frame(V(miserables)$name,miserables.clust)
colnames(df.clust) = c("V","clust")

#view graph colored by cluster
V(miserables)$label.color = miserables.clust
plot(miserables, vertex.shape="none")
#tkplot(miserables)

#=======================classification=======================================
source('classification.R')

class = classify.multiple(times = 100, graph = miserables, clusters = df.clust, distance = logforest_dist, alpha = 0.01)

#plot
V(miserables)$label.color = class
plot(miserables, vertex.shape="none")

#=====================modularity and error alpha-resistance==========================
alphas <- seq(0.01, 1.0, by=0.01)
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
        class = classify.multiple(times = 100, graph = miserables, clusters = df.clust, distance = dist.vect[[i]], alpha = a)
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
png(filename=file.path(folder,'miserables_mod2.png'))
plot(g.mods)
dev.off()

#view dependency of alpha and error
acc.melted.df <- melt(acc.df, id=c("alpha"))
g.acc<-ggplot(acc.melted.df) + 
    geom_point(aes(alpha, value, colour=variable)) +
    #geom_smooth(aes(alpha, value, colour=variable)) +
    geom_line(aes(alpha, value, colour=variable))
png(filename=file.path(folder,'miserables_error2.png'))
plot(g.acc)
dev.off()