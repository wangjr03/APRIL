#merge networks which share same TF
source("functions.R")
library(igraph)
library(Matrix)

args <- commandArgs(T)


n_cluster <- as.character(args[1])

if (n_cluster == "default"){
  n_cluster <- floor(dim(sub_graph_TF)[1]/20)
}else{
  n_cluster <- as.numeric(n_cluster)
}

#load subgraph
load("../output/network_adjacency_matrix.Rdata")

chrname <- paste0("chr",c(1:22,"X"))

#load island name
island <- read.table("../output/island_loc.txt")
island <- subset(island,island[,1]%in%chrname)
island[,1] <- factor(as.character(island[,1]),levels = chrname)

#load island-TF matrix
load("../output/frag_TF_mat.Rdata")

island_name <- apply(island,1,function(x){gsub(" ","",paste0(x,collapse = "-"))})

sub_graph_TF <- sapply(sub_adj_mat,function(x){
  
  tmp_name <- rownames(x)
  
  loc <- match(tmp_name,island_name)
  
  tmp_profile <- apply(frag_TF_mat[loc,],2,mean)
  return(tmp_profile)
  
})

sub_graph_TF <- t(sub_graph_TF)

n <- dim(sub_graph_TF)[1]

rownames(sub_graph_TF) <- 1:n

sub_sd_row <- apply(sub_graph_TF,1,sd)

sub_sd_col <- apply(sub_graph_TF,2,sd)

sub_graph_TF <- sub_graph_TF[which(sub_sd_row!=0),which(sub_sd_col!=0)]

sub_graph_TF <- apply(sub_graph_TF,2,function(x){x/sum(x)})

#cluster subgraph
get_skeleton <- function(adj,n_sub){
  
  n <- sum(n_sub)+1
  
  N <- dim(adj)[1]
  
  delete_id <- c()
  
  for(i in n:N){
    
    linked_nodes <- which(adj[i,]>0)
    
    id <- sapply(linked_nodes,function(x){min(which(cumsum(n_sub)>=x))})
    
    if(length(unique(id))==1){
      
      delete_id <- c(delete_id,i)
      
    }
  }
  if(length(delete_id)>0){
    
    return(adj[-delete_id,-delete_id])
    
  }else {
    
    return(adj)
    
  }
  
  
  
}

library(pheatmap)

sub_graph_cluster <- pheatmap(sub_graph_TF,clustering_distance_rows="correlation",clustering_distance_cols="correlation",clustering_method = "complete")

#generate new adj matrix


cluster <- cutree(sub_graph_cluster$tree_row,k=n_cluster )


merged_graph <- list()

m <- dim(frag_TF_mat)[2]

all_tf <- colnames(frag_TF_mat)

skeleton_graph <- list()

sub_adj_mat <- sapply(sub_adj_mat,function(x){1*as.matrix(x)})

for(iter in 1:n_cluster){
  
  component_id <- as.numeric(names(cluster)[which(cluster==iter)])
  
  component_adj <- sub_adj_mat[component_id]
  
  frag_id <- lapply(component_adj,rownames)
  
  tmp_frag_tf <- lapply(frag_id,function(x){
    
    loc <- match(x,island_name)
    
    frag_TF_mat[loc,] > 0
    
  })
  
  n <- sum(sapply(component_adj,dim)[1,])+m
  
  n_sub <- sapply(component_adj,dim)[1,]
  
  n_cum <- c(0,cumsum(n_sub))
  
  merged_adj <- matrix(0,ncol=n,nrow=n)
  
  node_name <- c(unlist(frag_id),all_tf)
  
  colnames(merged_adj) <- rownames(merged_adj) <- node_name
  
  N <- sum(sapply(component_adj,dim)[1,])
  
  for(i in 2:(length(n_sub)+1)){
    
    merged_adj[(n_cum[i-1]+1):n_cum[i],(n_cum[i-1]+1):n_cum[i]] <- component_adj[[i-1]]
    
    
    merged_adj[(N+1):n,(n_cum[i-1]+1):n_cum[i]] <- t(tmp_frag_tf[[i-1]])
    
    merged_adj[(n_cum[i-1]+1):n_cum[i],(N+1):n] <- tmp_frag_tf[[i-1]]
    
  }
  
  row_count <- apply(merged_adj,1,sum)
  
  merged_adj <- merged_adj[which(row_count>0),which(row_count>0)]
  
  diag(merged_adj) <- 0

  merged_graph[[iter]] <- merged_adj
  
  skeleton_graph[[iter]] <- get_skeleton(merged_adj,n_sub)
  
}

merged_graph <- c(merged_graph,sub_adj_mat[which(sub_sd_row==0)])

skeleton_graph <- c(skeleton_graph,sub_adj_mat[which(sub_sd_row==0)])

#a <- graph_from_adjacency_matrix(skeleton,mode="undirected")
#plot(a,layout=layout_with_lgl(a),vertex.size=8)

save(merged_graph,skeleton_graph,file="../output/TF_linked_subgraph.Rdata")


#annotate skeleton graph

load("../output/frag_annotation.Rdata")
sub_graph_nodes <- sapply(sub_adj_mat,rownames)
n<- length(sub_adj_mat)

island_anno <- do.call(rbind,lapply(1:n,function(x){data.frame(frag_name = sub_graph_nodes[[x]],sub_graph_id = x)}))
island_anno$anno <- frag_anno[match(island_anno[,1],island_name)]

n <- length(skeleton_graph)

skeleton_igraph <- list()

for(i in 1:n){
  
  tmp_graph <- skeleton_graph[[i]]
  
  node_name <- rownames(tmp_graph)
  
  sub_graph_id <- island_anno[match(node_name,island_anno[,1]),2]
  
  #use 0 for TF nodes
  sub_graph_id[which(is.na(sub_graph_id))] <- 0
  
  sub_graph_anno <- island_anno[match(node_name,island_anno[,1]),3]
  
  #use T for TF nodes
  sub_graph_anno[which(is.na(sub_graph_anno))] <- "T"
  
  tmp_igraph <- graph_from_adjacency_matrix(tmp_graph,mod="undirected")
  
  tmp_igraph <- set_vertex_attr(tmp_igraph, name='sub_graph_id', value=sub_graph_id)
  tmp_igraph <- set_vertex_attr(tmp_igraph, name='annotation', value=sub_graph_anno)
  
  skeleton_igraph[[i]] <- tmp_igraph
  
}

save(skeleton_igraph,n_cluster,file="../output/merged_skeleton_igraph.Rdata")
