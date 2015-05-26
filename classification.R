#function that marks one random vertex in each cluster
#and reclassifies according to distance matrix
classify <- function(D, marked, knn=1) {
    D.marked = D[,as.vector(marked$V)]
    
    graph.clust = sapply(rownames(D), function(v){
        min.dist = sort(D.marked[v,])[1:knn]
        min.names = names(min.dist)
        min.class = sapply(min.names, function(min.name) marked[marked$V==min.name,2])
        
        class = as.numeric(names(which(table(min.class)==max(table(min.class)))))
        if(length(class)>1) {
            df.min.class = data.frame(class = min.class, dist = min.dist)
            df.class.avgdist = aggregate(dist~class ,df.min.class[df.min.class$class %in% class,], mean)
            class = df.class.avgdist[which.min(df.class.avgdist$dist),]$class
        }
        class
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
        dim(class) <- c(dim(D)[1],i)
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

select.marked.multiple <- function(clusters, n=100, k=1) {
    lapply(1:n, function(i) {
        marked = as.vector(sapply(1:max(clusters$clust), function(i){
            sample(clusters[clusters$clust==i,1], k)
        }))
        df.marked = clusters[clusters$V %in% marked,]
    })
}



