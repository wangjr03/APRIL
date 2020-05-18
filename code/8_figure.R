
args <- commandArgs(T)

ct <- as.numerc(args[1])


#test for whole network, if the genes within same hub have higher correlation
load("../data/enh_act_mat.Rdata")
load("../data/gene_act_mat_updated.Rdata")

load("../output/island_annotation.Rdata")

gene_id <- sapply(gene_id,function(x){paste0(x,collapse=";")})

island_name <- data.frame(island_name,gene_id)

#gene's is more active

pdf("../figure/gene_exp_comparison.pdf",width = 3)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

col=c("#DB843D","#70588E")

l <- subset(gene_id,gene_id!="0")

l <- as.numeric(unlist(strsplit(l,";")))




boxplot(subset(gene_act_mat[l,ct],gene_act_mat[l,ct]>0),subset(gene_act_mat[-l,ct],gene_act_mat[-l,ct]>0),outline=F,border=col,names=NULL,ylab="Gene expression")
add_legend("topright", legend = c("Genes included in components","All expressed genes"), pch=20, 
           col=col,
           horiz=F, bty='n', cex=1)


dev.off()


g_g_corr <- list()

for(i in 1:length(skeleton_igraph)){
  
  tmp_g <- skeleton_igraph[[i]]
  
  p_id <- grep("P",V(tmp_g)$annotation)
  
  p_name <- V(tmp_g)$name[p_id]
  
  gene_id <- island_name[match(p_name,island_name[,1]),2]
  
  gene_id <- as.numeric(unlist(strsplit(as.character(gene_id),";")))
  
  gene_id <- gene_id[which(gene_act_mat[gene_id,ct]>0)]
  
  n <- length(gene_id)
  
  if(n<2){g_g_corr[[i]] <- NULL}else{
    
    corr <- cor(t(gene_act_mat[gene_id,]))
    
    corr <- corr[which(upper.tri(corr))]
    g_g_corr[[i]] <- corr
    
  }
  
  
}

n_gene <- as.numeric(unlist(strsplit(as.character(island_name[,2]),";")))




rd_corr <- function(x){
  
  rd_g <- sample(n_gene,2)
  
  corr <- cor(t(gene_act_mat[rd_g,]))
  
  corr <- corr[which(upper.tri(corr))]
  
  return(corr)
  
}

n_gene <- setdiff(as.numeric(unlist(strsplit(as.character(island_name[,2]),";"))),0)

n_gene <- n_gene[which(gene_act_mat[n_gene,ct]>0)]

generate_rd_corr <- function(x){
  
  rd_g <- sample(n_gene,2)
  
  corr <- cor(t(gene_act_mat[rd_g,]))
  
  corr <- corr[which(upper.tri(corr))]
  
  return(corr)
  
}

rd_corr <- sapply(1:1000,generate_rd_corr)


#bg corr: random pair across all genes

bg_index <- sapply(1:1000,function(x){sample(1:dim(gene_act_mat)[1],2)})


bg_corr <- apply(bg_index,2,function(x){cor(as.numeric(gene_act_mat[x[1],]),as.numeric(gene_act_mat[x[2],]))})


library("ggpubr")

names=c("Paired genes\n same components","Randomly paired\n connected genes","Randomly paired\nall genes")

w_score <- unlist(g_g_corr)

o_score <- unlist(rd_corr)

b_score <- na.omit(bg_corr)

df <- data.frame(Correlation=c(w_score,o_score,b_score),group=rep(names,c(length(w_score),length(o_score),length(b_score))))

df[,2] <- factor(as.character(df[,2]),levels = names)

#col=c("#DB843D","#4473A8","#70588E","#519BAD","#89A54D","#93A8CF","#AA4744")



pdf("../figure/Gene-gene_correlation.pdf",width = 5)


col=c("#FF0000","#DB843D","#70588E")

my_comparisons <- list(c(names[1],names[2]),c(names[1],names[3]))
ggboxplot(df, x = "group", y = "Correlation",col="group",palette=col)+ 
  stat_compare_means(comparisons = my_comparisons)

dev.off()


#repeat correlation analysis on enh-gene

load("../output/island_annotation.Rdata")

gene_id <- sapply(gene_id,function(x){paste0(x,collapse=";")})

enh_id <- sapply(enh_id,function(x){paste0(x,collapse=";")})


island_name <- data.frame(island_name,gene_id,enh_id)


e_g_corr <- list()

for(i in 1:length(skeleton_igraph)){
  
  tmp_g <- skeleton_igraph[[i]]
  
  p_id <- grep("P",V(tmp_g)$annotation)
  
  p_name <- V(tmp_g)$name[p_id]
  
  gene_id <- island_name[match(p_name,island_name[,1]),2]
  
  gene_id <- as.numeric(unlist(strsplit(as.character(gene_id),";")))
  
  gene_id <- gene_id[which(gene_act_mat[gene_id,ct]>0)]
  
  e_id <- grep("E",V(tmp_g)$annotation)
  
  e_name <- V(tmp_g)$name[e_id]
  
  enh_id <- island_name[match(e_name,island_name[,1]),3]
  
  enh_id <- as.numeric(unlist(strsplit(as.character(enh_id),";")))
  
  enh_id <- enh_id[which(enh_act_mat[enh_id,ct]>0)]
  
  cob <- expand.grid(enh_id,gene_id)
  
  n <- dim(cob)[1]
  
  if(n<2){e_g_corr[[i]] <- NULL}else{
    
    corr <- apply(cob,1,function(x){
      
      cor(as.numeric(gene_act_mat[as.numeric(x[2]),]),as.numeric(enh_act_mat[as.numeric(x[1]),]))
      
    })
    
    e_g_corr[[i]] <- corr
    
  }
  
  
}

n_gene <- setdiff(as.numeric(unlist(strsplit(as.character(island_name[,2]),";"))),0)

n_enh <- setdiff(as.numeric(unlist(strsplit(as.character(island_name[,3]),";"))),0)


rd_corr <- function(x){
  
  rd_g <- sample(n_gene,1)
  
  rd_e <- sample(n_enh,1)
  
  corr <- cor(as.numeric(gene_act_mat[rd_g,]),as.numeric(enh_act_mat[rd_e,]))
  
  
  return(corr)
  
}

rd_corr <- sapply(1:13419,rd_corr)


#bg corr: random pair across all genes

bg_index_g <- sapply(1:13419,function(x){sample(1:dim(gene_act_mat)[1],1)})

bg_index_e <- sapply(1:13419,function(x){sample(1:dim(enh_act_mat)[1],1)})

bg_index <- cbind(bg_index_g,bg_index_e)

bg_corr <- apply(bg_index,1,function(x){cor(as.numeric(gene_act_mat[x[1],]),as.numeric(enh_act_mat[x[2],]))})



library("ggpubr")

names=c("Paired EP\n same components","Randomly paired\n connected E & P","Randomly paired\nall E & P")

w_score <- unlist(e_g_corr)

o_score <- unlist(rd_corr)

b_score <- na.omit(bg_corr)

df <- data.frame(Correlation=c(w_score,o_score,b_score),group=rep(names,c(length(w_score),length(o_score),length(b_score))))

df[,2] <- factor(as.character(df[,2]),levels = names)

#col=c("#DB843D","#4473A8","#70588E","#519BAD","#89A54D","#93A8CF","#AA4744")



pdf("../figure/EP_correlation.pdf",width = 6)


col=c("#FF0000","#DB843D","#70588E")

my_comparisons <- list(c(names[1],names[2]),c(names[1],names[3]))
ggboxplot(df, x = "group", y = "Correlation",col="group",palette=col)+ 
  stat_compare_means(comparisons = my_comparisons)

dev.off()




pdf("../figure/enh_act_comp.pdf",width = 3)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

col=c("#DB843D","#70588E")

l <- subset(enh_id,enh_id!="0")

l <- as.numeric(unlist(strsplit(l,";")))




boxplot(subset(enh_act_mat[l,ct],enh_act_mat[l,ct]>0),subset(enh_act_mat[-l,ct],enh_act_mat[-l,ct]>0),outline=F,border=col,names=NULL,ylab="Enhancer activity")
add_legend("topright", legend = c("Enhancers included in components","All active Enhancers"), pch=20, 
           col=col,
           horiz=F, bty='n', cex=1)


dev.off()





#first get TF binding for each hub

load("../output/frag_TF_mat.Rdata")

N <- length(skeleton_igraph)

ppi <- read.table("../data/BIOGRID_PPI_TF_NAME.bed",sep="\t",colClasses = "character")

val_frac <- c()

for(i in 1:N){
  print(i)
  tmp_graph <- skeleton_igraph[[i]]
  
  frag_id <- na.omit(match(V(tmp_graph)$name,island_name))
  
  frag_tf <-apply(frag_TF_mat[frag_id,],1,function(x){names(which(x>0))})
  
  n <- sapply(frag_tf,length)
  
  frag_tf <- frag_tf[which(n>0)]
  
  #frag_tf <- names(frag_tf)[which(frag_tf>0)]
  
  if(length(frag_tf)>1){
    
    n <- length(frag_tf)
    
    id <- expand.grid(1:n,1:n)
    
    id <- subset(id,id[,1]!=id[,2])
    
    pot_pair <- apply(id,1,function(x){
      
      expand.grid(frag_tf[[as.numeric(x[1])]],frag_tf[[as.numeric(x[2])]])
      
    })
    
    pot_pair <- do.call(rbind,pot_pair)
    
    match_id <- prodlim::row.match(pot_pair,ppi)
    
    val_frac <- c(val_frac,1-mean(is.na(match_id)))}
  
}


#generate bg TF val frac

frag_name <- lapply(skeleton_igraph, function(x){
  
  tmp_name <- V(x)$name
  
  #tmp_gwas_hit <- V(x)$GWAS_hit
  
  #return(data.frame(tmp_name,tmp_gwas_hit))
  return(data.frame(tmp_name))
})

for(i in 1:length(frag_name)){
  
  frag_name[[i]]$index <- i
  
  
}

frag_name <- na.omit(do.call(rbind,frag_name))

frag_name <- frag_name[grep("chr",frag_name[,1]),]


rd_val_list <- list()
for(j in 1:1000){
  print(j)
  rd_frag <- frag_name
  
  rd_frag[,2] <- sample(rd_frag[,2])
  
  rd_frag <- split(rd_frag[,1],rd_frag[,2])
  
  
  rd_val_frac <- c()
  for(i in 1:N){
    #print(i)
    
    
    tmp_graph <- rd_frag[[i]]
    
    frag_id <- na.omit(match(tmp_graph,island_name))
    
    frag_tf <-apply(frag_TF_mat[frag_id,],1,function(x){names(which(x>0))})
    
    n <- sapply(frag_tf,length)
    
    frag_tf <- frag_tf[which(n>0)]
    
    #frag_tf <- names(frag_tf)[which(frag_tf>0)]
    
    if(length(frag_tf)>1){
      
      n <- length(frag_tf)
      
      id <- expand.grid(1:n,1:n)
      
      id <- subset(id,id[,1]!=id[,2])
      
      pot_pair <- apply(id,1,function(x){
        
        expand.grid(frag_tf[[as.numeric(x[1])]],frag_tf[[as.numeric(x[2])]])
        
      })
      
      pot_pair <- do.call(rbind,pot_pair)
      
      match_id <- prodlim::row.match(pot_pair,ppi)
      
      rd_val_frac <- c(rd_val_frac,1-mean(is.na(match_id)))}
    
  }
  
  rd_val_list[[j]] <- rd_val_frac
  
}

rd_mean_frac <- na.omit(sapply(rd_val_list,mean))




df <- data.frame(dist=c(mean(val_frac),mean(rd_mean_frac)),group=c("Components",rep("Random components",1)))

Mean <- c(mean(val_frac),mean(rd_mean_frac,na.rm=T))

SD <- c(0,sd(rd_mean_frac,na.rm=T))

col=c("#DB843D","#70588E")


pdf("../figure/ppi.pdf",width=4)

ggplot(df, aes(x = group, y = dist,fill=factor(group))) +
  geom_bar(stat = "identity",width = 0.8,cex=2)+
  ggtitle("PPI validated fraction") + ylab("Mean fraction") +
  # add 68% CI errorbar 
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2)+theme(axis.text=element_text(size=12,face="bold"),
                                                                            panel.background = element_rect(fill = "transparent"), # bg of the panel
                                                                            plot.background = element_rect(fill = "transparent", color = NA),
                                                                            axis.line = element_line(colour = "black")# bg of the plot
  )+scale_fill_manual(values = col)+theme(legend.position="top")

dev.off()


#check TF expression
tf <- sapply(skeleton_igraph,function(x){V(x)$name})
tf <- unlist(tf)
tf <- tf[-grep("chr",tf)]
tf_ensg <- as.character(subset(gene_annotation[,7],gene_annotation[,6]%in%tf))

tf_exp <- gene_act_mat[match(tf_ensg,promoter[,4]),ct]

bg_tf_exp <- gene_act_mat[-match(tf_ensg,promoter[,4]),ct]

w_score <- tf_exp

o_score <- bg_tf_exp

pdf("../figure/TF_exp.pdf",width = 3)
names=c("TF\nwithin the network","TF\noutside the network")


col=c("#DB843D","#70588E")

boxplot(w_score,o_score,border=col,outline=F,names = c("",""),cex.axis=0.7)



add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

col=c("#DB843D","#70588E")

#boxplot(subset(k,k>0),subset(gene_act_mat[-l,ct],gene_act_mat[-l,ct]>0),outline=F,border=col,names=NULL,ylab="Gene expression")
add_legend("topleft", legend = names, pch=20, 
           col=col,
           horiz=F, bty='n', cex=1)


dev.off()


#TF-gene correlation
gene_anno <- read.table("../data/gene_annotation_V19.txt")

promoter <- read.table("../data/promoter.txt")

load("frag_annotation.Rdata")
load("island_annotation.Rdata")
load("frag_TF_mat.Rdata")
load("../data/gene_act_mat_updated.Rdata")

pro <- grep("P",frag_anno)

n <- length(pro)


tg_corr <- c()

for(i in 1:n){
  
  tmp_ind <- pro[i]
  
  tmp_gene <- gene_id[[tmp_ind]]
  
  tmp_tf <- names(which(frag_TF_mat[tmp_ind,]>0))
  
  tmp_ensg <- as.character(subset(gene_anno[,7],gene_anno[,6]%in%tmp_tf))
  
  tmp_tf_id <- match(tmp_ensg,promoter[,4])
  
  if(any(!is.na(tmp_tf_id))&any(gene_act_mat[tmp_gene,ct]>1)){
    print(i)
    
    tmp_gene <- subset(tmp_gene,gene_act_mat[tmp_gene,ct]>1)
    
    tmp_tf_id <- na.omit(tmp_tf_id)
    
    tmp_pair <- expand.grid(tmp_tf_id,tmp_gene)
    
    tmp_cor <- apply(tmp_pair,1,function(x){
      
      cor(as.numeric(gene_act_mat[as.numeric(x[1]),]),as.numeric(gene_act_mat[as.numeric(x[2]),]))
      
    })
    
    tg_corr <- c(tg_corr,tmp_cor)
    
  }else{
    
    next
    
  }
  
}

all_tf <- colnames(frag_TF_mat)

rd_gene <- unlist(gene_id)[which(unlist(gene_id)>0)]

bg_pair <- expand.grid(all_tf_id,rd_gene)


n <- dim(bg_pair)[1]

bg_pair <- bg_pair[sample(1:n,1000),]

bg_corr <- apply(bg_pair,1,function(x){
  
  cor(as.numeric(gene_act_mat[as.numeric(x[1]),]),as.numeric(gene_act_mat[as.numeric(x[2]),]))
  
})



all_ensg <- as.character(subset(gene_anno[,7],gene_anno[,6]%in%all_tf))

all_tf_id <- match(all_ensg,promoter[,4])

act_gene <- which(gene_act_mat[,ct]>0)

all_pair <- expand.grid(all_tf_id,act_gene)

n <- dim(all_pair)[1]

sp_pair <- all_pair[sample(1:n,1000),]

sp_corr <- apply(sp_pair,1,function(x){
  
  cor(as.numeric(gene_act_mat[as.numeric(x[1]),]),as.numeric(gene_act_mat[as.numeric(x[2]),]))
  
})





library("ggpubr")

names=c("TF-gene\nwithin the network","TF-genenowithin the network","random pair")

w_score <- tg_corr

o_score <- bg_corr

i_score <- sp_corr

df <- data.frame(LogP=c(w_score,o_score,i_score),group=rep(names,c(length(w_score),length(o_score),length(i_score))))

df[,2] <- factor(as.character(df[,2]),levels = names)

#col=c("#DB843D","#4473A8","#70588E","#519BAD","#89A54D","#93A8CF","#AA4744")



pdf("../figure/TF_gene_exp_corr.pdf",width = 4)


col=c("#FF0000","#DB843D","#70588E")

ylim1 = boxplot.stats(df$LogP)$stats[c(1, 5)]

my_comparisons <- list(c(names[1],names[2]),c(names[1],names[3]))
ggboxplot(df, x = "group", y = "LogP",col="group",palette=col,outlier.shape = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+labs(y="-log(P)")

dev.off()


