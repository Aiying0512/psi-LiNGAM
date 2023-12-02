setwd("/Users/xsa5db/Documents/LiNGAM")
## Method Comparison
## PC; GES; ICA-LiNGAM; direct LiNGAM; psi-LiNGAM
library(R.matlab)
library(pcalg)
library(rlingam)

library(igraph)
library(bnlearn)
library(ppcor)
library(RBGL)

library(reticulate)
use_python("/usr/bin/python3")
lingam <- import("lingam")

df.empty <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df.empty) <- c("PC",'GES',"ICA-LiNGAM","direct-LiNGAM","psi-LiNGAM")
tprl <- df.empty
fprl <- df.empty
fdrl <- df.empty
SHDl <- df.empty
RFl <- c()
#RFl <- c()
df_empty <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_empty) <- c("p","PC",'GES',"ICA-LiNGAM","direct-LiNGAM","psi-LiNGAM")
mtpr <- df_empty
mfpr <- df_empty
mfdr <- df_empty
mSHD <- df_empty

shd_cal <- function(x,bn1){
  e <- empty.graph(as.character(1:p))
  amat(e,check.cycles=F) <- as.matrix(x)
  SHD <- bnlearn::shd(e,bn1)
  SHD
}



## number of variables i.e. # of nodes
##degree = 1,2,4
##p = 50,100,200 
n <- 500
plist <- c(50,100,200)
degree <- 2

N <- 5
h<- 2
p <- plist[h]

set.seed(1000)

for(t in 1:N){
# Let's make a simple simulation this time
rDAG <- randDAG(p,degree, method = "watts", DAG = TRUE)
#edge <- rDAG@edgeData@data
rgDAG <- graph_from_graphnel(rDAG)
### Get adjacency matrix
##row out
##column in
A <- as.matrix(get.adjacency(rgDAG))
num <- sum(A)
B <- A
B[A==1] <- runif(num, min = 0.3, max = 0.8)
#heatmap(B,Rowv = NA, Colv = NA)
bn1 <- empty.graph(as.character(1:p))
amat(bn1) <- A


###### Topologically sort the graph
###### Generate observations
top_ind <- as.numeric(tsort(rDAG))
X <- matrix(0,n,p)
for(i in 1:p){
  w <- B[,top_ind[i]]
  X[,top_ind[i]] <-  X%*%w + (rchisq(n,df=1)-1)/sqrt(2)
  #X[,top_ind[i]] <-  X%*%w + (rexp(n,rate =1)-1)
  
}


ts <- X
c <- p*(p-1)/2

### pc algorithm ########
gCPDAG <- pc(suffStat = list(C = cor(ts), n = n),
             indepTest = gaussCItest, ## (partial correlations)
             alpha = 0.05/c, p=p, verbose = FALSE)

adj_pc <- showAmat(gCPDAG)
adj_pc[adj_pc==1] <- 0
adj_pc[adj_pc!=0] <- 1
pcgp <- graph_from_adjacency_matrix(adj_pc)


## ges ####
score <- new("GaussL0penObsScore", ts)
ges.fit <- ges(score)

essg <- getGraph(ges.fit$essgraph)
g_ges <- graph_from_graphnel(essg)
adj_ges <- as.matrix(get.adjacency(g_ges))

## lingam ICA ####
lingam.fit <- lingam(ts)
g_lingam <- graph_from_adjacency_matrix(t(lingam.fit$Bpruned),mode = "directed",weighted = TRUE,diag = FALSE)
#plot.igraph(g_lingam,layout = mycoord,edge.color="orange")

adj_ica <- get.adjacency(g_lingam)

## direct lingam
model = lingam$DirectLiNGAM()
model$fit(ts)

adj_dl <- t(model$adjacency_matrix_)
adj_dl[adj_dl!=0] <-1

## psi lingam
#############################
####### psi-LiNGAM ##########
#############################

prior <- psi_prior(X)
sum(prior!=0&A==1)

model2 = lingam$DirectLiNGAM(prior_knowledge = prior)
model2$fit(ts)

adj_pdl <- t(model2$adjacency_matrix_)
adj_pdl[abs(adj_pdl)<0.2] <- 0
B1 <- adj_pdl
adj_pdl[adj_pdl!=0] <- 1


adj_est <- list(adj_pc, adj_ges, adj_ica, adj_dl, adj_pdl)

tpr_com <- lapply(adj_est, function(x) sum(A==1&x!=0)/sum(A)) %>% unlist
fpr_com <- lapply(adj_est, function(x) sum(A==0&x==0)/sum(A==0)) %>% unlist
fdr_com <- lapply(adj_est, function(x) 1-sum(A==1&x!=0)/sum(x!=0)) %>% unlist
shd_com <- lapply(adj_est, function(x) shd_cal(x,bn1)) %>% unlist
RF <- sqrt(sum((B - B1)*(B-B1)))/sqrt(sum(B1*B1))

RFl <- c(RFl,RF)

tprl[t,] <- tpr_com
fprl[t,] <- fpr_com
fdrl[t,] <- fdr_com
SHDl[t,] <- shd_com
}

#########
mtpr[h,] <- c(p,colMeans(tprl)) 
mfpr[h,] <- c(p,colMeans(fprl))
mfdr[h,] <- c(p,colMeans(fdrl))
mSHD[h,] <- c(p,colMeans(SHDl))

mean(RFl)
sd(RFl)
