###############################
###### psi-LiNGAM #############
###############################
psi_LiNGAM <- function(X,prior,mode=1,part_p=NULL){
  library(reticulate)
  use_python("/usr/bin/python3")
  lingam <- import("lingam")
  model2 = lingam$DirectLiNGAM(prior_knowledge = prior)
  model2$fit(X)
  if(mode==1){
    adj_pdl <- t(model2$adjacency_matrix_)
    adj_pdl[abs(adj_pdl)<0.2] <- 0
    B1 <- adj_pdl
    return(B1)
  }else{
  K <- model2$causal_order_ %>% unlist() +1
  
  Pr <- prior[K,K]
  Pr[Pr!=0] <-1
  Pr[lower.tri(Pr)]<-0
  pt <- abs(part_p)[K,K]
  pt[lower.tri(pt)] <- 0
  #Porder <- Pr[order(K),order(K)]
  #ts <- X
  ts <- X
  X1 <- ts[,K]
  B1 <- matrix(0,p,p)
  Pv <- matrix(0,p,p)
  p_adj <- matrix(0,p,p)
  
  nei <- ceiling(n/log(n))
  
  for(i in 1:p){
    if(sum(Pr[,i])>0){
      if(sum(Pr[,i])>(2*nei)){ind <- order(pt[,i],decreasing = T)[1:(nei*2)]}else{
        ind <- which(Pr[,i]==1)}
      lm.fit <- lm(X1[,i] ~ X1[,ind])
      B1[ind,i] <- summary(lm.fit)$coefficients[-1,1]
      Pv[ind,i] <- summary(lm.fit)$coefficients[-1,4]
      p_adj[ind,i] <- p.adjust(Pv[ind,i], n=p-1,method = "BH")
    }
  }
  
  #p.adjust(Pv[upper.tri(Pv)])
  
  B1 <- B1[order(K),order(K)]
  Pv <- Pv[order(K),order(K)]
  p_adj <- p_adj[order(K),order(K)]

  
  
  #adj_pdl[adj_pdl!=0] <- p.adjust(Pv[Pv!=0])
  
  B1[p_adj>0.05] <- 0
  return(list('B'=B1,'p_val'=p_adj))
  }
}
