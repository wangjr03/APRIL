######### effect size propagate #######
library(igraph)
library(stringr)

args <- commandArgs(T)
flag = as.character(args[1])

load("../output/merged_skeleton_igraph.Rdata")
load("../output/graph_effect.Rdata")

# remove the graph contains chrX
# remove_index <- sapply(skeleton_igraph, function(x){
#   tmp = grep(pattern = "chrX", x)
#   ifelse(length(tmp)==0, return(TRUE), return(FALSE))
# })
#
# skeleton_igraph = skeleton_igraph[remove_index]

##### step 1 #####
# label disease gene & controlled gene

# read gene annotation
promoter <- read.table("../data/promoter.txt", header = F)
promoter[,1] <- as.character(promoter[,1])
promoter[,4] <- as.character(promoter[,4])

gene_annotation <- read.table("../data/gene_annotation_V19.txt", header=T)
gene_annotation[,6] <- as.character(gene_annotation[,6])
gene_annotation[,7] <- as.character(gene_annotation[,7])

# map to the graph
get_gene <- function(x, y){ # x:graph, y:effect
  # get promoter index
  gene_index <- grep(pattern = "P", V(x)$annotation)
  
  # get gene name
  gene_pos <- V(x)$name[gene_index]
  tmp_split <- str_split(gene_pos, pattern = "-")
  gene_name <- lapply(tmp_split, function(y){
    tmp_start <- as.numeric(y[2])
    tmp_end <- as.numeric(y[3])
    tmp_pos <- promoter[promoter[,1] == y[1] & promoter[,2] < tmp_end & promoter[,3] > tmp_start,]
    tmp_index <- unique (grep(pattern = paste(tmp_pos[,4], collapse="|"), gene_annotation[,7]))
    return(gene_annotation[tmp_index,6])
  })
  
  # annotate disease_gene
  disease_label <- which( y[gene_index] == max(y[gene_index]) )  # highest rank
  
  record = list()
  record$gene_index <- gene_index
  record$gene_name<- gene_name
  record$disease_gene <- disease_label
  record$control_gene <- ifelse(length(gene_index)>1,
                                sample(c(1:length(gene_index))[-disease_label], 1),
                                0)
  return(record)
}

graph_gene <- mapply(get_gene, skeleton_igraph, graph_effect, SIMPLIFY = F)

# remove the graph contains less than 5 promoters
graph_index <- sapply(graph_gene, function(x) return( length(x$gene_index)>5) )
disease_graph <- skeleton_igraph[graph_index]
graph_disease_gene <- graph_gene[graph_index]
graph_disease_effect <- graph_effect[graph_index]


#### Step 2 ####
# HotNet to propagate effect, 
# HotNet algorithm: L_gamma^-1 %*% B
# define the L_gamma matrix, L = D - A; L_gamma = L + gamma*I

# calculate L_gamma matrix
L_matrix <- function(x){
  gam = mean(degree(x))
  tmp_adj <- as_adj(x)
  # get D
  tmp_D <- diag(apply(tmp_adj, 1, sum))
  tmp_L <- tmp_D - tmp_adj
  return(tmp_L + gam*diag(dim(tmp_L)[1]))
}

# HotNet: L_gamma^-1 %*% label
hotnet <- function(L,b,x){
  if (flag == "GENE") {
    b[ c(1:length(b))[-x$gene_index] ]=0
  }
  if (flag == "OTHER") {
    b[x$gene_index] = 0
  }
  
  final_label <- solve(L) %*% b 
  frag_name <- colnames(final_label)
  return(final_label)
}

# label propagation
L_gamma <- lapply(disease_graph, L_matrix)
result <- mapply(hotnet, L_gamma, graph_disease_effect, graph_disease_gene, SIMPLIFY = F)

# gene results
gene_result <- mapply(function(x,y) return(x[y$gene_index]),
                      result, graph_disease_gene)

#### step 3 ####
# draw performance boxplot and results

disease_rank <- mapply(function(x,y){
  tmp_order <- rank(x)
  tmp_order[y$disease_gene] / length(x)
}, gene_result, graph_disease_gene)


control_rank <- mapply(function(x,y){
  tmp_order <- rank(x)
  tmp_order[y$control_gene] / length(x)
}, gene_result, graph_disease_gene)

result_label <- c( rep("disease gene",length(disease_rank)), rep("controlled gene",length(control_rank)) )
result_rank <- c(disease_rank, control_rank)
result_dat <- data.frame(cbind(result_rank, result_label))
result_dat$result_rank <- as.numeric(as.character(result_dat$result_rank))
test_score <- pairwise.t.test(result_dat$result_rank, result_dat$result_label)
pdf("propagation_performance.pdf")
boxplot(result_rank~result_label, result_dat,border = c("#DB843D","#70588E"), outline=F, xlab = "",
        ylab = "Percentile rank", main = paste0("Pairwise t-test\n p value = ",test_score$p.value))
dev.off()

# rank in each network
pred_rank <- mapply(function(x,y) y$gene_name[x[order(x,decreasing = T)]], 
                    gene_result, graph_disease_gene, SIMPLIFY = F)

save(pred_rank, "../output/propagate_rank.Rdata")





