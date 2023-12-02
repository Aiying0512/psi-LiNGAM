#### Toy example
library(R.matlab)
library(igraph)
library(pheatmap)
library(reticulate)
use_python("/usr/bin/python3")
lingam <- import("lingam")


net <- readMat("sim1.mat")$net[1,,]
diag(net) <- 0

### set color and scale
library(RColorBrewer)
display.brewer.all()
color <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))(11)
library(scales)
show_col(color)

myBreaks <- c(seq(-0.2,0,0.1), 
              seq(0.1, 0.5, length.out=5))

myColor = color[-c(1:4)]

rownames(net) <- paste0('V',c(1:5))
colnames(net) <- paste0('V',c(1:5))

pheatmap(net,cluster_rows = F, cluster_cols = F, color  = myColor, breaks = myBreaks, 
         labels_row=paste0('V',c(1:5)),labels_col = paste0('V',c(1:5)), angle_col = 0,fontsize = 15 )
ts <- readMat("sim1.mat")$ts[1:10000,]
#writeMat("X1.mat",ts=ts)
p <- dim(ts)[2]
n <- dim(ts)[1]

g1 <- graph_from_adjacency_matrix(net,mode = "directed",weighted = TRUE,diag = FALSE)
#mycoord <- c(0,0,5,0,6.54,4.75,2.50,7.69,-1.54,4.75)
#mycoord<-matrix(mycoord,5,2,byrow=TRUE)
mycoord <- layout_in_circle(g1)
plot.igraph(g1,layout = mycoord, vertex.label=V(g1)$number,vertex.size=30,vertex.label.cex=2, 
            edge.arrow.size=1.5, edge.arrow.width=1.5,edge.width=4)

### pc algorithm ########
gCPDAG <- pc(suffStat = list(C = cor(ts), n = dim(ts)[1]),
             indepTest = gaussCItest, ## (partial correlations)
             alpha = 0.05, p=dim(ts)[2], verbose = FALSE)

adj_pc <- showAmat(gCPDAG)
adj_pc[adj_pc==1] <- 0
adj_pc[adj_pc!=0] <- 1
pcgp <- graph_from_adjacency_matrix(adj_pc)
plot.igraph(pcgp,layout = mycoord,edge.color="orange",vertex.size=30,vertex.label.cex=2, 
            edge.arrow.size=1.5, edge.arrow.width=1.5,edge.width=4)

## ges ####
score <- new("GaussL0penObsScore", ts)
ges.fit <- ges(score)

essg <- getGraph(ges.fit$essgraph)
g_ges <- graph_from_graphnel(essg)
adj_ges <- as.matrix(get.adjacency(g_ges))
plot.igraph(g_ges,layout = mycoord,edge.color="orange",vertex.size=30,vertex.label.cex=2, 
            edge.arrow.size=1.5, edge.arrow.width=1.5,edge.width=4)


## lingam ICA ####
lingam.fit <- lingam(ts)
g_lingam <- graph_from_adjacency_matrix(t(lingam.fit$Bpruned),mode = "directed",weighted = TRUE,diag = FALSE)
plot.igraph(g_lingam,layout = mycoord,edge.color="orange",vertex.size=30,vertex.label.cex=2, 
            edge.arrow.size=1.5, edge.arrow.width=1.5,edge.width=4)

adj_ica <- get.adjacency(g_lingam)

pheatmap(t(lingam.fit$Bpruned),cluster_rows = F, cluster_cols = F, color  = myColor, breaks = myBreaks, 
         labels_row=paste0('V',c(1:5)),labels_col = paste0('V',c(1:5)), angle_col = 0,fontsize = 15 )

### lingam direct ###
## direct lingam
model = lingam$DirectLiNGAM()
model$fit(ts)

adj_dl <- t(model$adjacency_matrix_)

g_dlingam <- graph_from_adjacency_matrix(adj_dl,mode = "directed",weighted = TRUE,diag = FALSE)
plot.igraph(g_dlingam,layout = mycoord,edge.color="orange",vertex.size=30,vertex.label.cex=2, 
            edge.arrow.size=1.5, edge.arrow.width=1.5,edge.width=4)

pheatmap(adj_dl ,cluster_rows = F, cluster_cols = F, color  = myColor, breaks = myBreaks, 
         labels_row=paste0('V',c(1:5)),labels_col = paste0('V',c(1:5)), angle_col = 0,fontsize = 15 )

### psi Lingam
source('psi_prior.R')
source('psi_lingam.R')

prior <- psi_prior(ts)$prior
adj_pdl <- psi_LiNGAM(ts,prior)

g_psilingam <- graph_from_adjacency_matrix(adj_pdl,mode = "directed",weighted = TRUE,diag = FALSE)
plot.igraph(g_psilingam,layout = mycoord,edge.color="orange",vertex.size=30,vertex.label.cex=2, 
            edge.arrow.size=1.5, edge.arrow.width=1.5,edge.width=4)

pheatmap(adj_pdl ,cluster_rows = F, cluster_cols = F, color  = myColor, breaks = myBreaks, 
         labels_row=paste0('V',c(1:5)),labels_col = paste0('V',c(1:5)), angle_col = 0,fontsize = 15 )
