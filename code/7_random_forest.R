############ use graph to find disease genes ##############
library(dplyr)
library(igraph)
library(stringr)
library(randomForest)
library(ROCR)
library(caret)

args <- commandArgs(T)
flag = as.numeric(args[1])
ntrees = as.numeric(args[2])
dis_gene = as.character(args[3])
sig_eqtl = as.character(args[4])

#### step1 ####
# get the TF features index for each nodes
# output: node_index: list, index in island_loc

# load graph, TF, annotation, fragments
load("../output/merged_skeleton_igraph.Rdata")
load("../output/frag_TF_mat.Rdata")
load("../output/island_annotation.Rdata")
frag <- read.table("../output/island_loc.txt")
frag[,1] <- as.character(frag[,1])

# remove the graph contains chrX
# remove_index <- sapply(skeleton_igraph, function(x){
#   tmp = grep(pattern = "chrX", x)
#   ifelse(length(tmp)==0, return(TRUE), return(FALSE))
# })
#
# skeleton_igraph = skeleton_igraph[remove_index]

# combine the name
frag_name <- sapply(1:nrow(frag), function(x){
  paste(frag[x,], collapse = "-")
})

# get the index of each fragments
node_index <- lapply(skeleton_igraph, function(x){
  match(V(x)$name, table = frag_name)
})


##### step 2 #######
# activity feature for enhancer and genes
load("../output/island_annotation.Rdata")

# enhancer activity features
enh <- read.table("../data/consensus_enhancer_coord.txt")
load("../data/enh_act_mat.Rdata")
# get activity score
enh_act <- enh_act_mat[,flag]
# calculate z-score
enh_mean <- apply(enh_act_mat, 1, mean)
enh_sd <- apply(enh_act_mat, 1, sd)
enh_z <- (enh_act - enh_mean) / enh_sd

# gene activity features
load("../data/gene_act_mat_updated.Rdata")
# get activity score
gene_act <- gene_act_mat[,flag]
# calculate z-score
gene_mean <- apply(gene_act_mat, 1, mean)
gene_sd <- apply(gene_act_mat, 1, sd)
gene_z <- (gene_act - gene_mean) / gene_sd


#### step3 ####
# annotate disease genes & negative genes
# output: graph_gene: list, index, gene_name, label(T or F)

# read gene annotation
promoter <- read.table("../data/promoter.txt", header = F)
promoter[,1] <- as.character(promoter[,1])
promoter[,4] <- as.character(promoter[,4])
gene_annotation <- read.table("../data/gene_annotation_V19.txt", header=T)
gene_annotation[,6] <- as.character(gene_annotation[,6])
gene_annotation[,7] <- as.character(gene_annotation[,7])

# read disease genes
disease_gene <- read.table(dis_gene)
disease_gene <- as.character(disease_gene[,1])
disease_gene <- gene_annotation[match(table = gene_annotation[,6], disease_gene),7]

# map to the graph
get_gene <- function(x,y){# x: graph, y:node_index
  # get promoter index
  gene_index <- grep(pattern = "P", V(x)$annotation)
  
  # get gene id
  gene_pos <- frag[y[gene_index],]
  gene_name <- apply(gene_pos, 1, function(z){
    tmp_start <- as.numeric(z[2])
    tmp_end <- as.numeric(z[3])
    tmp_pos <- promoter[promoter[,1] == z[1] & promoter[,2] < tmp_end & promoter[,3] > tmp_start,4][1]
    return(tmp_pos)
  })
  
  disease_label <- match(table = disease_gene, gene_name)
  disease_label <- ifelse( is.na(disease_label), 0, 1)
  
  # annotate neg gene, 0 is unlabelled
  if( sum(disease_label)>0 ){
    tmp_distance <- distances(x)
    sub_distance <- tmp_distance[gene_index[disease_label==1],gene_index] # row: disease, col: genes
    tmp_neg_select <- ifelse(is.matrix(sub_distance), which(apply(sub_distance,2,min) > 2), which(sub_distance > 2))
  }else{
    tmp_neg_select <- 1:length(disease_label)
  }
  disease_label[tmp_neg_select] <- -1
  
  # graph contain no gene
  if(length(gene_index)==0) disease_label <- NULL
  
  record = list()
  record$gene_index <- gene_index
  record$gene_id<- gene_name
  record$label <- disease_label
  return(record)
}

graph_gene <- mapply(get_gene, skeleton_igraph, node_index, SIMPLIFY = F)

load("../output/graph_effect.Rdata")

# eqtl set 
eqtl <- read.table(sig_eqtl)
eqtl[,3] <- as.character(eqtl[,3])

#### step 4 ####
# generate neg & pos features
data_feature <- function(x,y,z,w){ # x is subgraph, y is graph gene, z is node_index, w graph_effect
  # graph features
  tmp_degree <- degree(x)
  tmp_betweenness <- betweenness(x)
  tmp_closeness <- closeness(x)
  tmp_pagerank <- page.rank(x)$vector
  tmp_distance <- distances(x)
  gene_label <- y$label
  
  # get feature: several genes per graph
  tmp_feature <- mapply(function(k, p){ # multiple genes, gene_id
    # find neighbor
    neighbor_graph_index <- which(tmp_distance[k,]==1) # index in graph
    neighbor <- z[neighbor_graph_index] # index in frag_TF
    neighbor_graph_index <- neighbor_graph_index[!is.na(neighbor)] # remove TF nodes
    neighbor <- neighbor[!is.na(neighbor)] # remove TF nodes
    
    # neighbor effect size
    neighbor_effect_list <- w[neighbor_graph_index]
    neighbor_max_effect <- max(neighbor_effect_list)
    neighbor_sum_effect <- sum(neighbor_effect_list)
    n_neighbor <- length(neighbor)
    
    # average weight
    if (n_neighbor == 1){
      neighbor_feature <- frag_TF_mat[neighbor,]
    } else{
      if (neighbor_sum_effect == 0){
        neighbor_feature <- apply(frag_TF_mat[neighbor,], 2, mean)
      }else{
        neighbor_feature <- apply(frag_TF_mat[neighbor,], 2, function(TF_tmp){
          weighted.mean(TF_tmp, neighbor_effect_list)
        })
      }
    }

    ## eqtl
    tmp_eqtl <- eqtl[eqtl[,3]==p,]
    
    eqtl_p = 0.05
    eqtl_n = 0
    if(nrow(tmp_eqtl)>0){
      tmp_neighbor <- frag[neighbor,]
      for(i in 1:nrow(tmp_eqtl)){
        if(sum(tmp_neighbor[,1] == tmp_eqtl[i,1] & tmp_neighbor[,2] < tmp_eqtl[i,2] &
               tmp_neighbor[,3] > tmp_eqtl[i,2])>0){
          eqtl_p <- min(eqtl_p, tmp_eqtl[i,4])
          eqtl_n <- eqtl_n + 1
        }
      }
    }

    # their own features
    self_TF <- frag_TF_mat[z[k],]
    self_degree <- tmp_degree[k]
    self_betweenness <- tmp_betweenness[k]
    self_closeness <- tmp_closeness[k]
    self_pagerank <- tmp_pagerank[k]
    
    # gene z score and act
    self_act <- gene_act[gene_id[[z[k]]][1]]
    na_flag <- is.na(gene_z[gene_id[[z[k]]][1]])
    self_z <- ifelse(na_flag, 0, gene_z[gene_id[[z[k]]][1]]) 
    
    # neighbor z score and act
    enh_tmp_index <- enh_id[z[neighbor]]
    neighbor_act <- unlist(sapply(enh_tmp_index, function(x) enh_act[x]))
    neighbor_z <- unlist(sapply(enh_tmp_index, function(x) enh_z[x]))
    neighbor_act_stat <- summary(neighbor_act)
    neighbor_z_stat <- summary(neighbor_z)
    
    if(length(neighbor_act) == 0){ # no enhancer data
      neighbor_act_stat[] <- 0
      neighbor_z_stat[] <- 0
    }
    
    # feature rename
    names(neighbor_feature) <- paste0("neighbor_", names(neighbor_feature) )
    names(neighbor_max_effect) <- "max_neighbor_effect"
    names(neighbor_sum_effect) <- "sum_neighbor_effect"
    names(n_neighbor) <- "num_neighbor"
    names(self_TF) <- paste0("promoter_", names(self_TF) )
    names(self_degree) <- "degree"
    names(self_betweenness) <- "betweenness"
    names(self_closeness) <- "closeness"
    names(self_pagerank) <- "PageRank"
    names(neighbor_act_stat) <- paste0("neighbor_act_", names(neighbor_act_stat))
    names(neighbor_z_stat) <- paste0("neighbor_z_", names(neighbor_z_stat))
    names(self_act) <- "gene_act"
    names(self_z) <- "gene_z"
    names(eqtl_p) <- "eqtl_pvalue"
    names(eqtl_n) <- "eqtl_number"
    
    c(eqtl_n, eqtl_p, self_degree, self_betweenness, self_closeness, self_pagerank, neighbor_max_effect, neighbor_sum_effect,
      n_neighbor, self_act, self_z, neighbor_act_stat, neighbor_z_stat, self_TF, neighbor_feature)
  }, y$gene_index, y$gene_id, SIMPLIFY = TRUE)
  
  tmp_feature <- cbind(gene_label, t(tmp_feature))
  rownames(tmp_feature) <- y$gene_id
  return(tmp_feature)              
}

all_data <- mapply(data_feature, skeleton_igraph, graph_gene, node_index, graph_effect)
# remove null
null_index <- sapply(all_data, function(x) ncol(x) > 0 )
all_data <- all_data[null_index]
all_data <- do.call(rbind, all_data)

#### step4 ####
# train the model
pos <- all_data[all_data[,1]==1,-1]
neg <- all_data[all_data[,1]==-1,-1]

# balance data
if (nrow(neg) > nrow(pos)){
  neg <- neg[sample(1:nrow(neg), nrow(pos)),]
}else{
  pos <- pos[sample(1:nrow(pos), nrow(neg)),]
}

# label data
pos <- cbind(Label = 1, pos)
neg <- cbind(Label = 0, neg)
dat <- data.frame(rbind(neg,pos))
unlabel <- data.frame(all_data[all_data[,1]==0,])
dat$Label <- ifelse(dat$Label==1, "disease genes", "controlled genes")
dat$Label <- as.factor(dat$Label)

# cross-validation
folds = 10
val <- list(predictions = list(), labels = list(), importance = list())
dat_folds <- createFolds(1:nrow(dat), k = folds)
for (i in 1:folds) {
  test <- dat[dat_folds[[i]],]
  train <- dat[-dat_folds[[i]],]
  rf <- randomForest(Label~., train, ntree = ntrees)
  test_pred <- predict(rf, test[,-1], type = "prob")
  val$predictions <- append(val$predictions, list( as.vector(test_pred[,2]) ) ) # record probability
  val$labels <- append(val$labels, list( ifelse(test$Label == 'disease genes', 1, 0) ) ) # record true labels
  val$importance <- append(val$importance , list( as.vector(importance(rf)[,1]) ) ) # record importance
}

## Draw the averaged ROC curve and calculate AUC ##
pdf("../randomForest_diseaseGene.pdf")
ROCperf <- prediction(val$predictions, val$labels) %>% performance("tpr","fpr") # calculate TPR & FPR
ROCauc <- prediction(val$predictions, val$labels) %>% performance("auc")
(mean_auroc <- mean(unlist(ROCauc@y.values))) #calculate AUC
plot(ROCperf, col="grey",lty=6, main=paste("Mean AUC:",mean_auroc)) # individual ROC
plot(ROCperf, lwd=3, avg="vertical", add=T) # averaged ROC
dev.off()

## average feature importance
im <- matrix(unlist(val$importance), byrow = FALSE, ncol= folds)
im <- apply(im, 1, mean)     
names(im) <- colnames(dat)[-1]   
im <- data.frame(im[order(im, decreasing = T)])
write.table(im, "../output/importance_randomForest.txt", row.names = T, col.names = F, quote = F, sep = "\t")

## predict score            
prob_disease_test <- cbind(rownames(dat[unlist(dat_folds),]), unlist(val$prediction)) # test in 10-fold

rf <- randomForest(Label~., dat, ntree = ntrees)
test_pred <- predict(rf, unlabel[,-1], type = "prob")
prob_disease_unlabel <- test_pred[,2]
write.table(prob_disease, "../output/prob_disease_gene.txt", row.names = F, col.names = F, quote = F, sep = "\t")



