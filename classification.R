#function that marks one random vertex in each cluster
#and reclassifies according to distance matrix
classify <- function(D, clusters) {
    marked = sapply(1:max(clusters$clust), function(i){
        sample(clusters[clusters$clust==i,1], 1)
    })
    marked = levels(marked)[marked]
    df.marked = clusters[clusters$V %in% marked,]
    
    D.marked = D[,marked]
    
    graph.clust = sapply(rownames(D), function(v){
        marked.v = colnames(D.marked)[which.min(D.marked[v,])]
        df.marked[df.marked$V==marked.v,2]
    })
    graph.clust
}

classify.multiple <- function(times ,graph, clusters, distance, alpha) {
    D = distance(graph, alpha)
    
    #classification for numerous times
    class <- vector()
    for(i in seq(times)) {
        cat("=")
        class <- c(class,classify(D, clusters, distance))
        dim(class) <- c(77,i)
    }
    rownames(class) <- V(graph)$name
    
    #voting across all classifications
    class.vote <- apply(class,1,function(x) as.numeric(names(which.max(table(x)))))
    class.vote
}
