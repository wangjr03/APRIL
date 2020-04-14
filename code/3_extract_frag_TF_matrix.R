#for each fragment, calculate its TF profile based on enhancer/promoter
args <- commandArgs(T)

cell_index <- as.numeric(args[1])
TFthershold <- as.numeric(args[2])

#get enh_tf_matrix 
chr <- paste0("chr",c(1:22,"X"))

enh_coord <- read.table("../data/consensus_enhancer_coord.txt")

enh_tf_list <- list()

enh_use_pos <- list()

k <- 1

for(i in chr){
  
  load(paste0("../data/enh_tf_Rdata/",i,".Rdata"))
  
  enh_tf_list[[k]] <- enh_tf
  
  enh_use_pos[[k]] <- use_pos
  
  k <- k+1
}

#get promoter_tf_matrix 
chr <- paste0("chr",c(1:22,"X"))

promoter <- read.table("../data/promoter.txt")

promoter_tf_list <- list()

promoter_use_pos <- list()

k <- 1

for(i in chr){
  
  load(paste0("../data/promoter_tf_Rdata/",i,".Rdata"))
  
  promoter_tf_list[[k]] <- enh_tf
  
  promoter_use_pos[[k]] <- use_pos
  
  k <- k+1
}

#load subgraph

#load("../output/network_adjacency_matrix.Rdata")
island <- read.table("../output/island_loc.txt")
load("../output/island_annotation.Rdata")

#purify active enhancer and active gene

load("../data/enh_act_mat.Rdata")
load("../data/gene_act_mat_updated.Rdata")


act_enh_id <- sapply(enh_id,function(x){
  
  act_id <- subset(x,enh_act_mat[x,cell_index]>TFthershold)
  
  if(length(act_id)>0){return(act_id)}else{
    
    return(0)
    
  }
  
})

act_gene_id <- sapply(gene_id,function(x){
  
  act_id <- subset(x,gene_act_mat[x,cell_index]>TFthershold)
  
  if(length(act_id)>0){return(act_id)}else{
    
    return(0)
    
  }
  
})

act_gene_id <- gene_id

frag_anno <- c()

n <- dim(island)[1]

for(i in 1:n){
  
  if(!0%in%act_enh_id[[i]]& !0%in% act_gene_id[[i]]){
    
    frag_anno[i] <- "E/P"
    
    
  }else if(!0%in%act_enh_id[[i]]){
    
    frag_anno[i] <- "E"
    
  }else if(!0%in% act_gene_id[[i]]){
    
    frag_anno[i] <- "P"
    
  }else{
    
    frag_anno[i] <- "F"
    
  }
  
  
}

#generate frag_TF matrix
all_tf <- colnames(enh_tf)

n <- dim(island)[1] 

m <- length(all_tf)

frag_TF_mat <- matrix(0,ncol = m,nrow=n)

chr_order <- paste0("chr",c(1:22,"X"))

enh_coord[,1] <- factor(enh_coord[,1],levels=chr_order)

enh_chr_count <- as.numeric(table(enh_coord[,1]))

promoter[,1] <- factor(promoter[,1],levels=chr_order)

pro_chr_count <- as.numeric(table(promoter[,1]))

enh_chr_cum <- cumsum(c(0,enh_chr_count))

pro_chr_cum <- cumsum(c(0,pro_chr_count))

for(i in 1:n){
  
  tmp_chr <- which(chr==island[i,1])
  
  tmp_chr_id <- tmp_chr
  
  tmp_enh_id <- match(act_enh_id[[i]]-enh_chr_cum[tmp_chr_id],enh_use_pos[[tmp_chr]])
  
  tmp_gene_id <- match(act_gene_id[[i]]-pro_chr_cum[tmp_chr_id],promoter_use_pos[[tmp_chr]])
  
  tmp_profile <- rbind(enh_tf_list[[tmp_chr]][tmp_enh_id,],promoter_tf_list[[tmp_chr]][tmp_gene_id,])
  
  tmp_profile <- apply(tmp_profile,2,function(x){mean(x,na.rm = T)})
  
  tmp_profile[which(is.na(tmp_profile))] <- 0
  
  frag_TF_mat[i,] <- tmp_profile
}


#purify expressed TF
gene_name <- read.table("../data/gene_name.txt")

all_tf <- colnames(enh_tf)

gene_annotation <- read.table("../data/gene_annotation_V19.txt")

ensg_name <- as.character(gene_annotation[match(all_tf,gene_annotation[,6]),7])

expressed_tf <- which(gene_act_mat[match(ensg_name,gene_name[,1]),cell_index]>TFthershold)

frag_TF_mat <- frag_TF_mat[,expressed_tf]

colnames(frag_TF_mat) <- all_tf[expressed_tf]

save(frag_TF_mat,file="../output/frag_TF_mat.Rdata")

save(frag_anno,file="../output/frag_annotation.Rdata")
