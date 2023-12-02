##### Get Prior #######

#############################
####### psi-LiNGAM ##########
#############################
psi_prior <- function(X, thresh=NULL){
n <- dim(X)[1]
p <- dim(X)[2]

z_star2 <- NULL

if(n/p > 100){
  z_star2 <- qnorm(0.95)
}else{
  if(is.null(thresh)){thresh <- ifelse(n>=p, 0.4, 0.6)}
}
library(huge)
library(ppcor)
X <- huge.npn(X)
if(n>p){
  cor_partial <- pcor(X)
  if(is.null(z_star2)){
  z_star2 <- quantile(abs(cor_partial$statistic),thresh)}
  part_p <- cor_partial$estimate
  part_p[abs(cor_partial$statistic)<z_star2] <- 0
} else{
cor_pearson <- cor(X)
fisher_z <- 1/2 * sqrt(n-3)*log((1+cor_pearson)/(1-cor_pearson)) 
diag(fisher_z) <- 0

pho <- fisher_z
z_star1 <- quantile(abs(fisher_z),thresh)
pho[abs(pho)<z_star1] <- 0
#heatmap(pho,Rowv = NA, Colv = NA)

## latent transformation###
###choose neibourhood###
mid <- abs(pho)
nei <- ceiling(n/log(n))
max <- rep(0,p)
ind <- matrix(0,p,nei)
for(i in 1:p){
  max[i] <- min(sum(pho[i,]>0),nei)
  ind[i,] <-arrayInd(sort.list(mid[i,],decreasing=T)[1:nei],dim(mid))[,1]
}

###calculate precision matrix###
M <- data.frame(1:p,max)
#diag(pho) <-1
V <- matrix(0,p,p)
pz <- matrix(0,p,p)

###n<p run this one###
for (i in 1:(p-1)){
  for (j in (i+1):p){
    if(min(M[c(i,j),2])!=0){
      k <-c(i,j)[which.min(M[c(i,j),2])]
      
      idx <- ind[k,1:max[k]]
      if(length(which(idx==j))==0) {idx <- c(j,idx)} else {idx <- c(j,idx[-which(idx==j)])}
      if(length(which(idx==i))==0) {idx <- c(i,idx)} else {idx <- c(i,idx[-which(idx==i)])}
      
      
      sub <- pho[idx,idx]
      eigvals <- eigen(sub, only.values=T)$values
      perturb <- max(max(eigvals) - length(idx)*min(eigvals), 0)/(length(idx)-1)
      sub <- sub + diag(length(idx))*perturb
      inv <- solve(sub)
      diag <- diag(inv)
      V[i,j] <- inv[2]/sqrt(diag[1]*diag[2])
      V[j,i] <- V[i,j]
      #pz[i,j] <- (1+V[i,j])/(1-V[i,j])
      pz[i,j] <- sqrt(n-length(idx)-3)*0.5*log((1+V[i,j])/(1-V[i,j]))
      pz[j,i] <- pz[i,j]
      # phiz[i,j] <- qnorm(2*pnorm(pz[i,j])-1)
    }
  }
}


part_p <- pz
z_star2 <- quantile(abs(pz),thresh)
part_p[abs(part_p)<z_star2] <- 0
}


#sum(part_p&A==1)
#sum(pho&A==1)
#sum((part_p!=0|pho!=0)&A==1)

prior <- part_p
prior[prior!=0] <- -1

return(list('prior'=prior,'partial_cor' = part_p))
}


