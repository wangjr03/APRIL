#### effect size of nodes #####
library(igraph)
library(stringr)

args <- commandArgs(T)
SNP_file = as.character(args[1])

# read snp file
snp <- read.table(SNP_file, sep = "\t")
snp[,1] <- as.character(snp[,1])

# read fragment coordinate
frag <- read.table("../output/island_loc.txt")
frag[,1] <- as.character(frag[,1])

# match
result = c()
chr = paste0("chr", c(1:22,"X"))
for (name in chr) {
  frag_chr <- frag[frag[,1]==name,] # have sorted
  snp_chr <- snp[snp[,1] == name,]
  snp_chr <- snp_chr[order(snp_chr[,2]),]
  
  effect <- rep(0, nrow(frag_chr))
  
  frag_p = 1
  snp_p = 1
  frag_start <- frag_chr[frag_p,2]
  frag_end <- frag_chr[frag_p,3]
  
  while (snp_p <= nrow(snp_chr) & frag_p <= nrow(frag_chr)) {
    if(snp_chr[snp_p,2] < frag_start){
      snp_p <- snp_p + 1
      next
    }
    if(snp_chr[snp_p,2] >= frag_start & snp_chr[snp_p,2] <= frag_end){
      effect[frag_p] = max(effect[frag_p], snp_chr[snp_p,3])
      snp_p <- snp_p + 1
      next
    }
    if(snp_chr[snp_p,2] > frag_end){
      frag_p = frag_p + 1
      frag_start <- frag_chr[frag_p,2]
      frag_end <- frag_chr[frag_p,3]
      next
    }
  }
  
  frag_chr <- cbind(frag_chr, effect)
  result <- rbind(result, frag_chr)
}

frag_effect <- result

save(frag_effect, file = "../output/frag_effect.Rdata")

# annoate effect size of nodes in each network
load("../output/merged_skeleton_igraph.Rdata")

# remove the graph contains chrX
# remove_index <- sapply(skeleton_igraph, function(x){
#   tmp = grep(pattern = "chrX", x)
#   ifelse(length(tmp)==0, return(TRUE), return(FALSE))
# })
#
# skeleton_igraph = skeleton_igraph[remove_index]

get_effect <- function(x){
  # get fragments index
  node_tmp_index <- grep(pattern = "T", V(x)$annotation, invert = TRUE)
  TF_length <- length( grep(pattern = "T", V(x)$annotation) )
  
  # split node name
  node_pos <- V(x)$name[node_tmp_index]
  tmp_split <- str_split(node_pos, pattern = "-")
  
  # match effect size
  node_effect <- sapply(tmp_split, function(y){
    tmp_start <- as.numeric(y[2])
    tmp_pos <- frag_effect[frag_effect[,1] == y[1] & frag_effect[,2] == tmp_start,]
    return(tmp_pos[1,4])
  })
  
  node_effect <- c(node_effect,rep(0, TF_length))
  return(node_effect)
}

graph_effect <- lapply(skeleton_igraph, get_effect)

save(graph_effect, file = "../output/graph_effect.Rdata")
