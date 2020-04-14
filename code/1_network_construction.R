#1 parameter:
#1: input 3D link data in 5 column format: chr frag1.start frag1.end frag2.start frag2.end
args <- commandArgs(T)

links <- read.table(args[1])

links <- subset(links,links[,1]%in%paste0("chr",c(1:22,"X")))


n <- dim(links)[1]

links <- cbind(links[,c(1,2,3)],links[,c(1,4,5)])

links$link_id <- 1:n

pt1 <- links[,c(1,2,3,7)]

pt2 <- links[,c(4,5,6,7)]

colnames(pt2) <- colnames(pt1)

links <- rbind(pt1,pt2)


#find hubs
#sort file based on chr and start first
links <- links[order(links[,1],links[,2]),]
n <- dim(links)[1]
pointer <- 1

tmp_hub <- 1

hub_id <- c()

island <- list()

island_id <- 1

pointer_frag <- links[1,]

tmp_island_id <- 1

i <- 1

hub_id <- 1

tmp_frag <- links[i,]

tmp_link_id <- tmp_frag[,4]

loc <- which(links[,4]==tmp_link_id)

loc <- subset(loc,loc>i)

hub_id[loc] <- tmp_hub

print("overlapping fragments...")

while(i<n){
  
  #move to next fragment
  i <- i+1
  
  tmp_frag <- links[i,]
  
  #check overlapping
  
  #if not on the same chr, no overlapping, save previous island location, island id+1, hub_id +1
  if(tmp_frag[,1]!=pointer_frag[,1]){
    
    island[[tmp_island_id]] <- pointer_frag
    
    tmp_island <- links[i,]
    
    tmp_island_id <- tmp_island_id + 1
    
    island_id[i] <- tmp_island_id
    
    tmp_hub <- tmp_hub + 1
    
    
    #if already asigned hub_id, then don't touch it
    if(is.na(hub_id[i])){
    hub_id[i] <- tmp_hub
    
    tmp_link_id <- tmp_frag[,4]
    
    loc <- which(links[,4]==tmp_link_id)
    
    loc <- subset(loc,loc>i)
    
    hub_id[loc] <- tmp_hub
    
    }
    
    pointer_frag <- tmp_frag
    
  }else{
    
    #check overlapping if on the same chr
    
    if(tmp_frag[,2]<pointer_frag[,3]){
      
      #archive hub and island info for i th frag
      
      if(is.na(hub_id[i])){
      hub_id[i] <- tmp_hub}
      
      island_id[i] <- tmp_island_id
      
      #update hub_id for another side
      
      tmp_link_id <- links[i,4]
      
      loc <- which(links[,4]==tmp_link_id)
      
      loc <- subset(loc,loc>i)
      
      if(length(loc)>0){
        
        if(is.na(hub_id[loc])){
        
      hub_id[loc] <- tmp_hub}}
      
      
      #update pointer
      pointer_frag[,3] <- max(pointer_frag[,3],tmp_frag[,3])

      island[[tmp_island_id]] <- pointer_frag
    }else{
      
      #no overlapping
      island[[tmp_island_id]] <- pointer_frag
      
      tmp_island <- links[i,]
      
      tmp_island_id <- tmp_island_id + 1
      
      island_id[i] <- tmp_island_id
      
      tmp_hub <- tmp_hub + 1
      
      
      #if already asigned hub_id, then don't touch it
      if(is.na(hub_id[i])){
        hub_id[i] <- tmp_hub
        tmp_link_id <- tmp_frag[,4]
        
        loc <- which(links[,4]==tmp_link_id)
        
        loc <- subset(loc,loc>i)
        
        hub_id[loc] <- tmp_hub
        
      }
      
      pointer_frag <- tmp_frag
      
    
    }
  }
  
}

island[[tmp_island_id]] <- pointer_frag

links <- cbind(links,island_id,hub_id)

island <- do.call(rbind,island)[,c(1,2,3)]

island_name <- apply(island,1,function(x){gsub(" ","",paste0(x,collapse = "-"))})

#construct adjacancy matrix

n <- range(links$island_id)[2]

print("generating adjacency matrix")

require(Matrix)

adj_mat <- sparseMatrix(i=1,j=1,dims=c(n,n))

adj_mat[1,1] <- FALSE

for(i in 1:n){
  
  tmp_link_id <- subset(links$link_id,links$island_id==i)
  
  alt_loc <- links$island_id[which(links$link_id%in%tmp_link_id)]
  
  adj_mat[i,alt_loc] <- TRUE
  
}

rownames(adj_mat) <- colnames(adj_mat) <- island_name

##### detect hubs from it
source("functions.R")

sub_graph <- all_subgraph(adj_mat)

sub_adj_mat <- lapply(sub_graph,function(x){
  
  loc <- match(x,colnames(adj_mat))
  
  return(adj_mat[loc,loc])
  
})



save(sub_adj_mat,file="../output/network_adjacency_matrix.Rdata")
save(sub_graph,file="../output/network_nodes_name.Rdata")
write.table(island,"../output/island_loc.txt",col.names = F,row.names = F,sep="\t",quote=F)

