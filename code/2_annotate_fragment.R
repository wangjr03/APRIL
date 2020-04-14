#annotate fragment with enhancers or promoter or both
island <- read.table("../output/island_loc.txt")

enh_coord <- read.table("../data/consensus_enhancer_coord.txt")

promoter <- read.table("../data/promoter.txt")

#promoter <- promoter[,c(2,3,4,1)]

chr <- paste0("chr",c(1:22,"X"))


enh_coord[,1] <- factor(as.character(enh_coord[,1]),levels = chr)
promoter[,1] <- factor(as.character(promoter[,1]),levels = chr)
island[,1] <- factor(as.character(island[,1]),levels = chr)

#island <- island[order(island[,1],island[,2]),]
promoter <- promoter[order(promoter[,1],promoter[,2]),]
enh_coord <- enh_coord[order(enh_coord[,1],enh_coord[,2]),]

enh_coord[,1] <- as.character(enh_coord[,1])
promoter[,1] <- as.character(promoter[,1])
island[,1] <- as.character(island[,1])

# 
# get_match <- function(island,annotation){
#   
#   anno_pointer <- 1
#   
#   n <- dim(annotation)[1]
#   
#   island_pointer <- 1
#   
#   id <- c()
#   
#   m <- dim(island)[1]
#   
#   while(anno_pointer < n ){
#     
#     tmp_island <- island[island_pointer,]
#     
#     tmp_anno <- annotation[anno_pointer,]
#     
#     if(tmp_anno[,1]>tmp_island[,1]){
#       
#       id[island_pointer] <- 0
#       
#       island_pointer <- island_pointer+1
#       
#       if(island_pointer>m) break
#       
#     }else{
#       
#       if(tmp_anno[,1]<tmp_island[,1]){
#         
#         anno_pointer <- anno_pointer + 1
#         
#       }else{
#         
#         if(tmp_anno[,3] < tmp_island[,2]){
#           
#           anno_pointer <- anno_pointer + 1
#           
#         }else{
#           
#           if(tmp_anno[,3] > tmp_island[,2]){
#             
#             
#             id[island_pointer] <- 0
#             
#             island_pointer <- island_pointer+1
#             
#             if(island_pointer>m) break
#             
#           }else{
#             
#             id[island_pointer] <- anno_pointer
#             
#             island_pointer <- island_pointer + 1
#             
#             if(island_pointer>m) break
#             
#           }
#         }
#       }
#     }
#   }
#   
#   return(id)
# }

get_match <- function(x,y){
  
  tmp_chr <- as.character(x[1])
  
  tmp_start <- as.numeric(x[2])
  
  tmp_stop <- as.numeric(x[3])
  
  id <- which(y[,1]==tmp_chr&y[,2]<tmp_stop&y[,3]>tmp_start)
  
  if(length(id)==0){return(0)}else{
    
    return(id)
    
  }
  
}

enh_id <- apply(island,1,function(x){get_match(x,enh_coord)})
gene_id <- apply(island,1,function(x){get_match(x,promoter)})

save(enh_id,gene_id,file="../output/island_annotation.Rdata")

