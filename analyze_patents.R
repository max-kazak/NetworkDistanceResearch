library(igraph)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#prepare for parallel computation
library(doParallel)
library(foreach)
numWorkers=3

source('metrics.R')
folder = as.character(Sys.Date())
dir.create(folder, showWarnings = FALSE)

#Load Les Miserables graph
load("data//patents/patents_net_63-95_GB.Rdata")

#add fixed layout
mylayout <- layout.kamada.kawai(G.patents)
G.patents$layout <- mylayout

#choose palette
palette(brewer.pal(7, "Dark2"))
#display.brewer.all()

#view graph
plot(G.patents, vertex.label=NA, vertex.size=2)

#======================clusterization=======================================
#clusters in df
df.clust = data.frame(V(G.patents)$name,v.patents.cat)
colnames(df.clust) = c("V","clust")

#view graph colored by cluster
V(G.patents)$color = v.patents.cat

png(filename=file.path(folder,'patents_graph_clust.png'),
    width     = 3000,
    height    = 3000)
plot(G.patents, vertex.label=NA, vertex.size=2)
dev.off()

#=======================classification=======================================
source('classification.R')
mark.k = 10
knn.k = 3

l.marked.mult = select.marked.multiple(df.clust, k = mark.k)

strt <- Sys.time()
class = classify.multiple.voted(graph = G.patents, 
                                marked.mult = l.marked.mult, 
                                distance = logforest_dist, 
                                alpha = 0.01,
                                knn=knn.k)
print(Sys.time()-strt)

#plot
V(G.patents)$color = class
png(filename=file.path(folder,paste('patents_graph_class_m=',mark.k,'_k=',knn.k,'.png',sep='')),
    width     = 3000,
    height    = 3000)
plot(G.patents, vertex.label=NA, vertex.size=2)
dev.off()

#=====================modularity and error alpha-resistance==========================
alphas <- c(
    seq(0.001, 0.009, by=0.004), 
    seq(0.01,0.09, by=0.02),
    seq(0.1, 0.9, by=0.1)
    ,seq(0.91, 0.99, by=0.02)
    ,seq(0.991, 0.999, by=0.004) 
)
dist.vect <- c(plainwalk_dist, walk_dist , plainforest_dist, logforest_dist, communicability_dist, logcommunicability_dist)
names(dist.vect) <- c("plain_walk", "walk", "plain_forest", "log_forest", "communicability", "log_communicability")

cl = makeCluster(numWorkers)
clusterExport(cl, c("G.patents", "alphas", "l.marked.mult", "knn.k",
                    "classify.multiple", "classify",
                    "getT","spectral_radius", "addnames", "general_dist","plainwalk_dist", "walk_dist" , "plainforest_dist", 
                    "logforest_dist", "communicability_dist", "logcommunicability_dist"))
registerDoParallel(cl)

#calculate class for each dist*alpha*multiple
strt <- Sys.time()
l.dist.alpha.class <- lapply(dist.vect, function(dist) {
    print("processing")
    #calculate class for each alpha in alphas
    l.alpha.class = foreach(i = 1:length(alphas), .packages="igraph") %dopar% {
        cat("#")
        classify.multiple(G.patents, l.marked.mult, dist, alphas[i], knn = knn.k)        
    }     
})
print(Sys.time()-strt)
stopCluster(cl)

#calculate accuracy and modularity for each dist*alpha*multiple
l.dist.alpha.accmod <- lapply(l.dist.alpha.class, function(l.alpha.class) {
    l <- lapply(l.alpha.class, function(m.class) 
        apply(m.class,2, function(class) c(modularity(G.patents, class), 
                                           sum(class==v.patents.cat)/length(class)))
    )
    names(l) <- alphas
    l
})

#calculate quantiles of accuracy and modularity for each dist*alpha
l.dist.alpha.accmod.quant <- lapply(l.dist.alpha.accmod, function(l.alpha.accmod)
    lapply(l.alpha.accmod, function(m.accmod) {
        m.accmod.quant <- t(apply(m.accmod,1, function(accmod) quantile(accmod, probs=c(.25, .5, .75))))
        df.accmod.quant <- data.frame(m.accmod.quant)
        names(df.accmod.quant)<-c('q25','q50','q75')
        df.accmod.quant[["var"]]<-c('mod','acc')
        df.accmod.quant
    })
)

df.dist.alpha.accmod.quant<-melt(l.dist.alpha.accmod.quant, id=c('q25','q50','q75'))[,c(1,2,3,5,6,7)]
names(df.dist.alpha.accmod.quant) <- c('q25','q50','q75','measure','alpha','dist')
df.dist.alpha.accmod.quant$measure <- as.factor(df.dist.alpha.accmod.quant$measure)
df.dist.alpha.accmod.quant$dist <- as.factor(df.dist.alpha.accmod.quant$dist)


#==============visualize df.dist.alpha.accmod.quant============
ydim_mod = c(round(min(df.dist.alpha.accmod.quant[df.dist.alpha.accmod.quant$measure=='mod',]$q25),digits=1),
             round(max(df.dist.alpha.accmod.quant[df.dist.alpha.accmod.quant$measure=='mod',]$q75),digits=1))
ydim_acc = c(round(min(df.dist.alpha.accmod.quant[df.dist.alpha.accmod.quant$measure=='acc',]$q25),digits=1),
             round(max(df.dist.alpha.accmod.quant[df.dist.alpha.accmod.quant$measure=='acc',]$q75),digits=1))

#view dependency of alpha and modularity
g.mods <- ggplot(df.dist.alpha.accmod.quant[df.dist.alpha.accmod.quant$measure=='mod',], aes(x = alpha, group = dist)) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill=dist), alpha = .25) +
    geom_line(aes(y = q50, colour=dist), size=1) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Modularity") + ylab("modularoty") +
    scale_y_continuous(breaks = c(seq(from = ydim_mod[1], to = ydim_mod[2], by = 0.01)),
                       labels = c(seq(from = ydim_mod[1], to = ydim_mod[2], by = 0.01)))

ggsave(
    paste('patents_mod_m=',mark.k,'_k=',knn.k,'.png',sep=''),
    g.mods,
    path=folder,
    width = 10,
    height = 10
)

#view dependency of alpha and error
g.acc <- ggplot(df.dist.alpha.accmod.quant[df.dist.alpha.accmod.quant$measure=='acc',], aes(x = alpha, group = dist)) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill=dist), alpha = .25) +
    geom_line(aes(y = q50, colour=dist), size=1) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Error rate") + ylab("accuracy")+
    scale_y_continuous(breaks = c(seq(from = ydim_acc[1], to = ydim_acc[2], by = 0.01)),
                       labels = c(seq(from = ydim_acc[1], to = ydim_acc[2], by = 0.01)))

ggsave(
    paste('patents_acc_m=',mark.k,'_k=',knn.k,'.png',sep=''),
    g.acc,
    path=folder,
    width = 10,
    height = 10
)




