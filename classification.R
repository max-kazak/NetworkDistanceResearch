#function that marks one random vertex in each cluster
#and reclassifies according to distance matrix
classify <- function(graph, distance, alpha) {
    marked = sapply(1:max(df.clust$clust), function(i){
        sample(df.clust[df.clust$clust==i,1], 1)
    })
    marked = levels(marked)[marked]
    df.marked = df.clust[df.clust$V %in% marked,]
    
    D = distance(graph, alpha)
    D.marked = D[,marked]
    
    graph.clust = sapply(rownames(D), function(v){
        marked.v = colnames(D.marked)[which.min(D.marked[v,])]
        df.marked[df.marked$V==marked.v,2]
    })
    graph.clust
}

classify.multiple <- function(times ,graph, distance, alpha) {
    #classification for numerous times
    class <- vector()
    for(i in seq(times)) {
        class <- c(class,classify(graph, distance, alpha))
        dim(class) <- c(77,i)
    }
    rownames(class) <- V(graph)$name
    
    #voting across all classifications
    class.vote <- apply(class,1,function(x) as.numeric(names(which.max(table(x)))))
    class.vote
}