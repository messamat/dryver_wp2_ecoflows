#library(parallel)
#library(betapart)
library(dplyr)

divLeinster <- function(spxp, Z=NULL, q=2, check = TRUE){
  #Calcul the diversity of each site of sites by species matrix.
  #spxp columns and Z rows and columns are assumed to be in the same order.
  if (is.null(Z)) Z <- diag(ncol(spxp))
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  Zp <- Z %*% t(spxp)
  
  if (q != 1 & q != Inf){
    mat <- t(spxp) * (Zp)^(q-1)
    mat[is.na(mat)] <- 0
    D <- colSums(mat) ^ (1/(1-q))
  }
  if (q==Inf)  {
    D <- 1/ apply(Zp, 2, max)
  }
  if (q == 1){
    D <- apply(Zp^t(spxp), 2, function(x) 1/prod(x))
  }
  return(D)
}


abgDecompQ <- function(spxp, Z=diag(ncol(spxp)), q=2, site.weight = NULL, check=TRUE) {
  #Calcul the diversity of each site of sites by species matrix.
  #spxp columns and Z rows/cols are assumed to be in the same order.
  
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" do not have matching dimensions")}
    if (nrow(spxp) != length(site.weight) & !is.null(site.weight)){
      stop("object \"site.weight\" and object \"spxp\" do not have matching dimensions")
    }
  }
  if (is.null(site.weight)){
    site.weight <- rep(1/nrow(spxp), nrow(spxp))
  }else{
    site.weight <- site.weight/sum(site.weight)
  }
  
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  
  gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"))
  
  Gamma <- divLeinster(t(as.matrix(gamma.ab)), Z=Z , q=q, check = FALSE)
  Alphas <- divLeinster(spxp, Z=Z , q=q, check = FALSE)
  
  if (q != 1 & q != Inf) {
    mAlpha <- (sum(site.weight * (Alphas ^ (1 - q))))^(1 / (1 - q))
  }
  if (q==1){
    mAlpha <- exp(sum(site.weight * log(Alphas)))
  }
  if (q==Inf){
    mAlpha <- min(Alphas)
  }
  Beta <- Gamma / mAlpha
  
  names(Alphas) <- row.names(spxp)
  res <- list(Gamma=Gamma, Beta=Beta, mAlpha=mAlpha, Alphas=Alphas)
  
  return(res)
}

betapart.core.abund <- function(x){
  
  # test for a numeric matrix or data.frame
  if(! is.matrix(x)){
    x<-as.matrix(x)
  }
  
  if(! is.numeric(x))
    stop("The data in x is not numeric.",call.=TRUE)
  
  # simple test for positive data
  xvals <-  unique(as.vector(x))
  if (any(xvals<0)) 
    stop("The table contains negative values: data should be abundances.", call. = TRUE)
  
  pair.shared.abund<-matrix(nrow=nrow(x),ncol=nrow(x))
  rownames(pair.shared.abund)<-rownames(x)
  colnames(pair.shared.abund)<-rownames(x)
  
  pair.not.shared.abund<-matrix(nrow=nrow(x),ncol=nrow(x))
  rownames(pair.not.shared.abund)<-rownames(x)
  colnames(pair.not.shared.abund)<-rownames(x)
  
  for(i in 1:nrow(x)) {
    for(j in i:nrow(x)) {
      pair.shared.abund[j,i]<-sum(pmin(x[i,],x[j,]))
      pair.not.shared.abund[i,j]<-sum(x[i,])-sum(pmin(x[i,],x[j,]))
      pair.not.shared.abund[j,i]<-sum(x[j,])-sum(pmin(x[i,],x[j,]))
    }
  }
  
  pair.shared.abund<-as.dist(pair.shared.abund)
  pair.max.not.shared.abund<-pmax(as.dist(t(upper.tri(pair.not.shared.abund)*pair.not.shared.abund)), as.dist(pair.not.shared.abund))
  pair.min.not.shared.abund<-pmin(as.dist(t(upper.tri(pair.not.shared.abund)*pair.not.shared.abund)), as.dist(pair.not.shared.abund))
  
  multiple.shared.abund<-sum(x)-sum(apply(x,2,max))
  
  
  
  computations<-list(data=x, multiple.shared.abund=multiple.shared.abund, pair.shared.abund=pair.shared.abund, pair.min.not.shared.abund=pair.min.not.shared.abund, 
                     pair.max.not.shared.abund=pair.max.not.shared.abund, pair.sum.not.shared.abund=pair.min.not.shared.abund+pair.max.not.shared.abund)
  class(computations)<-"betapart.abund"
  
  return(computations)
} 


beta.multi.abund<-function (x, index.family = "bray") {
  index.family <- match.arg(index.family, c("bray", "ruzicka"))
  if (!inherits(x, "betapart.abund")) {
    x <- betapart.core.abund(x)
  }
  maxbibj <- sum(x$pair.max.not.shared.abund)
  minbibj <- sum(x$pair.min.not.shared.abund)
  switch(index.family, bray = {
    beta.bray.bal <- minbibj/(minbibj + x$multiple.shared.abund)
    beta.bray.gra <- (x$multiple.shared.abund/(minbibj + x$multiple.shared.abund)) * ((maxbibj - minbibj)/((2 * 
                                                                                                              x$multiple.shared.abund) + maxbibj + minbibj))
    beta.bray <- (minbibj + maxbibj)/(minbibj + maxbibj + 
                                        (2 * x$multiple.shared.abund))
    multi <- list(beta.BRAY.BAL = beta.bray.bal, beta.BRAY.GRA = beta.bray.gra, 
                  beta.BRAY = beta.bray)
  }, ruzicka = {
    beta.ruz.bal <- (2 * minbibj)/((2 * minbibj) + x$multiple.shared.abund)
    beta.ruz.gra <- (x$multiple.shared.abund/((2 * minbibj) + x$multiple.shared.abund)) * ((maxbibj - 
                                                                                              minbibj)/((x$multiple.shared.abund) + maxbibj + minbibj))
    beta.ruz <- (minbibj + maxbibj)/(minbibj + maxbibj + 
                                       x$multiple.shared.abund)
    multi <- list(beta.RUZ.BAL = beta.ruz.bal, beta.RUZ.GRA = beta.ruz.gra, 
                  beta.RUZ = beta.ruz)
  })
  return(multi)
}

# Alpha diversity of all patches
# alpha_div_raw <- function(metacom){
#   alpha_temp <- matrix(NA, ncol = length(metacom), nrow = nrow(metacom[[1]]))
#   for (tmp in 1:length(metacom)){
#     alpha_temp[,tmp] <- vegan::diversity(t(metacom[[tmp]]), MARGIN = 2, index = "invsimpson")
#   }
#   alpha_temp[is.infinite(alpha_temp)] <- 0
#   return(alpha_temp)
# }

# Alpha of wet communities and dry patches still occupied 
alpha_div <- function(metacom, flow_int){
  alpha_temp <- matrix(NA, ncol = length(metacom), nrow = nrow(metacom[[1]]))
  for (tmp in 1:length(metacom)){
    d_patches <- colnames(subset(flow_int, select = c(flow_int[tmp,]==0)))
    alpha <- vegan::diversity(t(metacom[[tmp]]), MARGIN = 2, index = "invsimpson")
    names(alpha) <- as.numeric(colnames(flow_int))
    alpha[names(alpha) %in% d_patches & is.infinite(alpha) == T] <- NA
    alpha_temp[,tmp] <- alpha
  }
  
  alpha_temp[is.infinite(alpha_temp)] <- 0
  return(alpha_temp)
}

# Alpha of wet communities and dry patches still occupied 
# alpha_div_unbal <- function(metacom, flow_int){
#   alpha_temp <- matrix(NA, ncol = length(metacom), nrow = nrow(metacom[[1]]))
#   for (tmp in 1:length(metacom)){
#     d_patches <- colnames(subset(flow_int, select = c(flow_int[tmp,]==0)))
#     alpha <- vegan::diversity(t(metacom[[tmp]]), MARGIN = 2, index = "invsimpson")
#     names(alpha) <- as.numeric(colnames(flow_int))
#     alpha[names(alpha) %in% d_patches & is.infinite(alpha) == T] <- NA
#     alpha_temp[,tmp] <- alpha
#   }
#   
#   alpha_temp[is.infinite(alpha_temp)] <- 0
#   return(alpha_temp)
# }


spatial_to_temporal <- function(metacom){
  metacom_temp <- lapply(1:nrow(metacom[[1]]), matrix, data =NA, nrow=length(metacom), ncol=ncol(metacom[[1]]))
  for(tmp in 1:length(metacom)){
    for(site in 1:nrow(metacom[[1]])){
      metacom_temp[[site]][tmp,] <- metacom[[tmp]][site,]
    }
  }
  return(metacom_temp)
}

# Weighting metacom by patch area
metacom_wg <- function(metacom, patch_size){
  patch_size <- matrix(rep(patch_size, ncol(metacom)), ncol = ncol(metacom), 
                       nrow = length(patch_size), byrow=F)
  metacom_wg <- metacom/patch_size
}


# beta_div <- function(metacom, patch_size){
#   metacom <- lapply(metacom, metacom_wg)
#   nb_cores <- detectCores() - 1 
#   cl <- makeCluster(nb_cores)
#   t = Sys.time()
#   beta <- parSapply(cl, metacom, FUN = beta.multi.abund)
#   stopCluster(cl)
#   Sys.time()-t
# }

# Beta div (Chao 2016)
# multi.bray.chao <- function(x){
#   
#   # test for a numeric matrix or data.frame
#   if(! is.matrix(x)){
#     x<-as.matrix(x)
#   }
#   
#   if(! is.numeric(x))
#     stop("The data in x is not numeric.",call.=TRUE)
#   
#   # simple test for positive data
#   xvals <-  unique(as.vector(x))
#   if (any(xvals<0)) 
#     stop("The table contains negative values: data should be abundances.", call. = TRUE)
#   
#   species.means<-apply(x,2,sum)/nrow(x)
#   
#   difs<-abs(x-t(matrix(rep(species.means,nrow(x)), nrow=ncol(x))))
#   
#   delta.BCnorm<-sum(difs)/(2*(1-1/nrow(x))*sum(x))
#   return(delta.BCnorm)
# }

# For moving averages
ma = function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}

# corr beta dissimil et dyst hydro
# weightd_met <- lapply(metacom, metacom_wg,  patch_size = hydro$patch_size$area)
# 
# d_dist <- as.dist(dendr_dist)
# cor_beta_dist <- function(metacom, dendr_dist){
#     metacom_dis <- vegan::vegdist(metacom, "bray")
#     rsq <- cor(metacom_dis, as.dist(dendr_dist), method = "spearman")
# }
# 
# rsq_met <- do.call(rbind, lapply(weightd_met, cor_beta_dist, dendr_dist = d_dist))
# hist(rsq_met)
# plot(rsq_met)


Anodiv<-function(spxp, structure, phy = NULL, weight = NULL, check = T) {
  #spxp : matrice with species as columns and sites as lines.
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}  
    if (!is.null(phy) & !inherits(phy, "phylo")){
      stop("object \"phy\" is not a phylo object")}
  }
  spxp <- sweep(spxp, 1, rowSums(spxp), "/") 
  
  if (!is.null(phy)){
    chao.objs <- chaoObjects(spxp, phy)
    spxp <- chao.objs$pi
    Z <- chao.objs$Z
  }else{
    Z <- diag(1, ncol(spxp))
  }
  #Check of the structures
  if (dim(structure)[1] != nrow(spxp) | dim(structure)[2] != 2){stop("error : 'structure' doesn't have the right dimensions")}
  
  if (is.null(weight)) {
    site.weight <- rep(1/nrow(spxp), nrow(spxp))
  } else {
    site.weight <- weight / sum(weight)
  }
  gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"))
  
  # gamma diversity #sp.weight
  Gamma <- divLeinster2(t(as.matrix(gamma.ab)), Z = Z , q = 1, check = FALSE)
  
  # Mean Alpha diversity
  Alphas <- divLeinster(spxp, Z = Z , q = 1, check = FALSE) 
  mAlpha <- exp(sum(site.weight * log(Alphas)))
  
  #Mean Alphas of intermediate classifications
  mAlphas.int <- apply(structure, 2, function(x){
    spxp.int<-aggregate(sweep(spxp, 1, site.weight, "*"), by=list(x), sum)
    site.weight.int <- aggregate(site.weight, by=list(x), sum)[,2]
    spxp.int <- sweep(spxp.int[,-1], 1, site.weight.int, "/")
    
    Alphas.int<- divLeinster(spxp.int, Z = Z, q = 1, check = FALSE)
    mAlpha.int <- exp(sum(site.weight.int * log(Alphas.int)))
  })
  names(mAlphas.int) <- c()  
  # list of results
  res<-c(Gamma = Gamma,
         Beta1 = Gamma/mAlphas.int[1],
         Beta2 = Gamma/mAlphas.int[2],
         Beta1x2 = prod(mAlphas.int)/(Gamma*mAlpha),
         mAlpha = mAlpha)
  return(res)
}

divLeinster2 <- function(spxp, Z=NULL, q=2, check = TRUE){
  #Calcul the diversity of each site of sites by species matrix. 
  #spxp columns and Z rows and columns are assumed to be in the right order.
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}  
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}  
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }  
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  Zp <- Z %*% t(spxp)
  
  if (q != 1 & q != Inf){
    mat <- t(spxp) * (Zp)^(q-1)
    mat[is.na(mat)] <- 0
    D <- colSums(mat) ^ (1/(1-q))
  }
  if (q == Inf)  {
    D <- 1/ apply(Zp, 2, max)
  }  
  if (q == 1){
    D <- apply(Zp^t(spxp), 2, function(x) 1/prod(x))
  }
  return(D)
}

chaoObjects <- function(spxp, phy){
  require(ape)
  require(phangorn)
  #Spxp has site in lines and species in columns. Within is relative abundances. 
  #Species are in the same order as tip.label of the phylogeny.
  
  Ancestors.sp <- lapply(1:length(phy$tip.label), function(x) c(Ancestors(phy, x, "all"), x))
  Branches <- lapply(Ancestors.sp, function(x) which((phy$edge[,1] %in% x) & (phy$edge[,2] %in% x)))
  Li <- unlist(lapply(Branches, function(x) sum(phy$edge.length[x]))) #Tip - root distances
  
  ultra <- all.equal.numeric(var(unlist(Li)), 0, tolerance = 1e-7) #Is it ultrametric
  if (!ultra) stop ("tree must be ultrametric.")
  
  freq.dummy <- unlist(lapply(1:length(Branches),function(i) rep(i,length(Branches[[i]]))))
  desc <- lapply(unlist(Branches), function(x) Descendants(phy, phy$edge[x,2], type ="tips")[[1]])
  
  tmp <- lapply(desc, function(desc_i){
    x <- rep(0, length(freq.dummy))
    x[freq.dummy %in% desc_i] <- 1
    return(x)
  })
  Z <- do.call(rbind, tmp)
  pi <- sweep(spxp[,freq.dummy], 2, phy$edge.length[unlist(Branches)], FUN = "*")/Li[1]
  
  res <- list(Z = Z, 
              pi = pi)
  res
}

# HierAnodiv <- function(spxp, structure, phy = NULL, weight = NULL, check = T, q = 1) {   # Before 22/11
#   #require(germinal)
#   #spxp : matrice with species as columns and sites as lines.
#   if (check){
#     if (!inherits(spxp, "matrix")) {
#       stop("object \"spxp\" is not of class \"matrix\"")}
#     # if (!inherits(Z, "matrix")) {
#     #   stop("object \"Z\" is not of class \"matrix\"")}
#     # if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
#     #   stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
#     if (!is.null(phy) & !inherits(phy, "phylo")){
#       stop("object \"phy\" is not a phylo object")}
#   }
#   spxp <- sweep(spxp, 1, rowSums(spxp), "/")
#   
#   if (!is.null(phy)){
#     chao.objs <- chaoObjects(spxp, phy)
#     spxp <- chao.objs$pi
#     Z <- chao.objs$Z
#   }else{
#     Z <- diag(1, ncol(spxp))
#   }
#   #Check of the structures
#   if (dim(structure)[1] != nrow(spxp) | dim(structure)[2] != 2){stop("error : 'structure' doesn't have the right dimensions")}
#   
#   if (is.null(weight)) {
#     site.weight <- rep(1/nrow(spxp), nrow(spxp))
#   } else {
#     site.weight <- weight / sum(weight)
#   }
#   gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"))
#   
#   # gamma diversity #sp.weight
#   Gamma <- divLeinster(t(as.matrix(gamma.ab)), Z = Z , q = q, check = FALSE)
#   
#   # Mean Alpha diversity
#   Alphas <- divLeinster(spxp, Z = Z , q = q, check = FALSE)
#   mAlpha <- exp(sum(site.weight * log(Alphas)))
#   if (q != 1 & q != Inf) {
#     mAlpha <- (sum(site.weight * (Alphas^(1 - q))))^(1/(1 -
#                                                           q))
#   }
#   if (q == 1) {
#     mAlpha <- exp(sum(site.weight * log(Alphas)))
#   }
#   if (q == Inf) {
#     mAlpha <- min(Alphas)
#   }
#   # Abg by site
#   tab_by_site <- do.call(rbind, lapply(unique(structure[,1]), function(i){
#     c(unlist(abgDecompQ(spxp[structure[,1] == i,], q= q)[1:3]), nsite = sum(structure[,1] == i), name = i)
#   }))
#   
#   spxp.int<-aggregate(sweep(spxp, 1, site.weight, "*"), by=list(structure[,1]), sum)
#   site.weight.int <- aggregate(site.weight, by=list(structure[,1]), sum)[,2]
#   spxp.int <- sweep(spxp.int[,-1], 1, site.weight.int, "/")
#   
#   Gamma.int<- divLeinster(spxp.int, Z = Z, q = 1, check = FALSE)
#   if (q != 1 & q != Inf) {
#     meanGamma.int <- (sum(site.weight * (Gamma.int^(1 - q))))^(1/(1 -q))
#   }
#   if (q == 1) {
#     meanGamma.int <- exp(sum(site.weight * log(Gamma.int)))
#   }
#   if (q == Inf) {
#     meanGamma.int <- min(Gamma.int)
#   }
#   # list of results
#   res<-list(c(Gamma = Gamma,
#               Beta1 = Gamma/meanGamma.int,
#               Beta2in1 = meanGamma.int/mAlpha,
#               mAlpha = mAlpha), tab_by_site =tab_by_site )
#   return(res)
# }

HierAnodiv <- function(spxp, structure, phy = NULL, weight = NULL, check = T, q = 1) {
  #require(germinal)
  #spxp : matrice with species as columns and sites as lines.
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    # if (!inherits(Z, "matrix")) {
    #   stop("object \"Z\" is not of class \"matrix\"")}
    # if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
    #   stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
    if (!is.null(phy) & !inherits(phy, "phylo")){
      stop("object \"phy\" is not a phylo object")}
  }
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  
  if (!is.null(phy)){
    chao.objs <- chaoObjects(spxp, phy)
    spxp <- chao.objs$pi
    Z <- chao.objs$Z
  }else{
    Z <- diag(1, ncol(spxp))
  }
  #Check of the structures
  if (dim(structure)[1] != nrow(spxp) | dim(structure)[2] != 2){stop("error : 'structure' doesn't have the right dimensions")}
  
  if (is.null(weight)) {
    site.weight <- rep(1/nrow(spxp), nrow(spxp))
  } else {
    site.weight <- weight / sum(weight)
  }
  gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"))
  
  # gamma diversity #sp.weight
  Gamma <- divLeinster(t(as.matrix(gamma.ab)), Z = Z , q = q, check = FALSE)
  
  # Mean Alpha diversity
  Alphas <- divLeinster(spxp, Z = Z , q = q, check = FALSE)
  
  if (q != 1 & q != Inf) {
    mAlpha <- (sum(site.weight * (Alphas^(1 - q))))^(1/(1 -
                                                          q))
  }
  if (q == 1) {
    mAlpha <- exp(sum(site.weight * log(Alphas)))
  }
  if (q == Inf) {
    mAlpha <- min(Alphas)
  }
  # Abg by site
  tab_by_site <- do.call(rbind, lapply(unique(structure[,1]), function(i){
    c(unlist(abgDecompQ(spxp[structure[,1] == i,,drop=FALSE], q= q)[1:3]), 
      nsite = sum(structure[,1] == i), name = i)
  }))
  
  spxp.int<-aggregate(sweep(spxp, 1, site.weight, "*"), by=list(structure[,1]), sum)
  site.weight.int <- aggregate(site.weight, by=list(structure[,1]), sum)[,2]
  spxp.int <- sweep(spxp.int[,-1], 1, site.weight.int, "/")
  
  Gamma.int<- divLeinster(spxp.int, Z = Z, q = 1, check = FALSE)
  if (q != 1 & q != Inf) {
    meanGamma.int <- (sum(site.weight.int * (Gamma.int^(1 - q))))^(1/(1 -q))
  }
  if (q == 1) {
    meanGamma.int <- exp(sum(site.weight.int * log(Gamma.int)))
  }
  if (q == Inf) {
    meanGamma.int <- min(Gamma.int)
  }
  # list of results
  res<-list(c(Gamma = Gamma,
              Beta1 = Gamma/meanGamma.int,
              Beta2in1 = meanGamma.int/mAlpha,
              mAlpha = mAlpha), tab_by_site =tab_by_site )
  return(res)
}
