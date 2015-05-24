#function that marks one random vertex in each cluster
#and reclassifies according to distance matrix
classify <- function(D, marked) {
    D.marked = D[,as.vector(marked$V)]
    
    graph.clust = sapply(rownames(D), function(v){
        marked.v = colnames(D.marked)[which.min(D.marked[v,])]
        marked[marked$V==marked.v,2]
    })
    graph.clust
}

classify.multiple <- function(graph, marked.mult, distance, alpha) {
    D = distance(graph, alpha)
    
    #classification for numerous times
    class <- vector()
    for(i in seq(length(marked.mult))) {
        #cat("=")
        class <- c(class,classify(D, marked.mult[[i]]))
        dim(class) <- c(77,i)
    }
    rownames(class) <- V(graph)$name
    
    class
}

classify.multiple.voted <- function(graph, marked.mult, distance, alpha) {
    class <- classify.multiple(graph, marked.mult, distance, alpha)
    #voting across all classifications
    class.vote <- apply(class,1,function(x) as.numeric(names(which.max(table(x)))))
    class.vote
}

select.marked.multiple <- function(clusters, num=100) {
    l.marked = list()
    for(i in seq(num)) {
        marked = sapply(1:max(clusters$clust), function(i){
            sample(clusters[clusters$clust==i,1], 1)
        })
        marked = levels(marked)[marked]
        df.marked = clusters[clusters$V %in% marked,]
        l.marked[[i]] = df.marked
    }
    l.marked
}



