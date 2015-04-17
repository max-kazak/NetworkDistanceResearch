

general_dist <- function(H_0, log = F) {
    n=dim(H_0)[1]
    Ones = matrix(1, ncol = n, nrow = n)
    
    if(log==T) H = log(H_0)
    if(log==F) H = H_0
    
    
    h = diag(diag(H))
    
    D = (h %*% t(Ones) + Ones %*% t(h) - H - t(H))/2
    
    return(D)
}

walk_dist <- function(G, alpha) {
    
    A = get.adjacency(G)
    t = getT(alpha, A)
    n = dim(A)[1]
    I = diag(n)
    
    H_0 =  solve(I - t*A)
    
    D = general_dist(H_0, log = T)
    
    return(D)
}


logforest_dist <- function(G ,alpha) {
    A = get.adjacency(G)
    L = graph.laplacian(G)
    t = getT(alpha, A)
    n = dim(A)[1]
    I = diag(n)
    
    H_0 =  solve(I + t*L)
    
    D = general_dist(H_0, log = T)
    
    return(D)    
}

plainforest_dist <- function(G, alpha) {
    A = get.adjacency(G)
    L = graph.laplacian(G)
    t = getT(alpha, A)
    n = dim(A)[1]
    I = diag(n)
    
    H =  solve(I + t*L)
    
    D = general_dist(H)
    
    return(D)
}

plainwalk_dist <- function(G, alpha) {
    A = get.adjacency(G)
    t = getT(alpha, A)
    n = dim(A)[1]
    I = diag(n)
    
    H =  solve(I - t*A)
    
    D = general_dist(H)
    
    return(D)
}

communicability_dist <- function(G, alpha) {
    A = get.adjacency(G)
    t = getT(alpha, A)
    
    H =  as.matrix(exp(t*A))
    
    D = general_dist(H)
    
    return(D)
}

logcommunicability_dist <- function(G, alpha) {
    A = get.adjacency(G)
    t = getT(alpha, A)
    
    H_0 =  as.matrix(exp(t*A))
    
    D = general_dist(H_0, log=T)
    
    return(D)
}


getT <- function(alpha, A) {
    ro = spectral_radius(A)
    t = 1/(1/alpha+ro)
    
    return(t)
}

spectral_radius <- function(M) return(max(abs(eigen(M)$values)))