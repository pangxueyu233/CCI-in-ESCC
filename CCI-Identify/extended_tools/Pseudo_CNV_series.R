Pseudo_CNV_OBJ <- function(Seurat_obj=Seurat_obj,sel_cells=epi_cells,ctrl_cells=ctrl_cells,
	sel_gene_pos=sel_gene_pos,anno_pos_need=anno_pos_need,project=project,save_path=save_path,file_prefix=file_prefix){
  require(Seurat)
  require(future.apply)
  require(trqwe)
  all_counts <- as.data.frame(GetAssayData(Seurat_obj,slot="counts"))
  sel_data <- all_counts[,sel_cells]
  ctrl_data <- all_counts[,ctrl_cells]
  ctrl_data_mean <- as.data.frame(rowSums(ctrl_data)/ncol(ctrl_data))
  colnames(ctrl_data_mean) <- "Standard_exp"
  epi_nor_data <- future_apply(sel_data,2, function(x) (x+1)/(ctrl_data_mean$Standard_exp+1))

  ##50 genes a bin
  all_chr_50b <- function(i){
    tmp <- subset(sel_gene_pos,chr_anno==i)
    tmp$order <- 1:nrow(tmp)
    rownames(tmp) <- tmp$ensembl_gene_id
    seed <- c(0,30,60,90)
    seed_names <- c("1st","2ed","3rd","4th")
    all_genes_plus_add_all <- c()
    for (order in 1:length(seed)){
      sel_order <- seed[order]
      if (nrow(tmp) > sel_order){
        if (sel_order != 0){
          sel_tmp <- tmp[1:sel_order,]
          sel_tmp$new_CNV <- paste("CNV",seed_names[order],"50bins_gap",sep="_")
          others_tmp <- tmp[setdiff(rownames(tmp),rownames(sel_tmp)),]
          }else{
            others_tmp <- tmp
            sel_tmp <- c()
        }
        }else{
          others_tmp <- tmp
          sel_tmp <- c()
      }
      if (nrow(others_tmp) > 50){
        tmp_order <- c()
        all_genes <- c()
        for(num in 1:floor(nrow(others_tmp)/50)){
          num_genes <- num*50
          order_genes <- rownames(others_tmp)[1:num_genes]
          order_genes <- setdiff(order_genes,tmp_order)
          tmp_genes <- others_tmp[order_genes,]
          tmp_genes$new_CNV <- paste("CNV",seed_names[order],"50bins",num,sep="_")
          tmp_order <- union(tmp_order,rownames(tmp_genes))
          all_genes <- rbind(all_genes,tmp_genes)
        }
        rownames(others_tmp) <- others_tmp$ensembl_gene_id
        if (length(setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id))==0){
          tmp_remind <- c()
          }else{
          tmp_remind <- others_tmp[setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id),]
          tmp_remind$new_CNV <- paste("CNV",seed_names[order],"50bins",round(nrow(others_tmp)/50),sep="_")
          }
        all_genes_plus <- rbind(all_genes,tmp_remind)
        all_genes_plus <- rbind(all_genes_plus,sel_tmp)
        }else{
          all_genes_plus <- others_tmp
          all_genes_plus$new_CNV <- paste("CNV",seed_names[order],"50bins",round(nrow(tmp)/50),sep="_")
          all_genes_plus <- rbind(all_genes_plus,sel_tmp)
      }
      rownames(all_genes_plus) <- all_genes_plus$ensembl_gene_id
      add_new <- data.frame(all_genes_plus[,ncol(all_genes_plus)],row.names=rownames(all_genes_plus))
      colnames(add_new) <- paste("CNV",seed_names[order],"50bins",sep="_")
      all_genes_plus <- cbind(all_genes_plus,add_new)
      all_genes_plus_sel <- as.data.frame(all_genes_plus[,setdiff(colnames(all_genes_plus),rownames(all_genes_plus_add_all))])
      if (ncol(all_genes_plus_sel)==1){
        colnames(all_genes_plus_sel) <- paste("CNV",seed_names[order],"50bins",sep="_")
        }else{
          all_genes_plus_sel <- all_genes_plus_sel
      }
      all_genes_plus_sel <- as.data.frame(t(all_genes_plus_sel))
      colnames(all_genes_plus_sel) <- rownames(all_genes_plus)
      all_genes_plus_add_all <- rbind(all_genes_plus_add_all,all_genes_plus_sel)
      #message(i,seed_names[order], " is done")
    }
    all_genes_plus_all_chr <- as.data.frame(t(all_genes_plus_add_all))
    return(all_genes_plus_all_chr)
    message(i," is done")
}

  all_chr_merge_50b <- future_lapply(unique(sel_gene_pos$chr_anno),all_chr_50b)
  all_genes_plus_all_chr_merge_50b <- do.call(rbind,all_chr_merge_50b)
  all_genes_plus_all_chr_merge_50b$CNV_1st_50 <- paste(all_genes_plus_all_chr_merge_50b$CNV_1st_50bins,all_genes_plus_all_chr_merge_50b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_50b$CNV_2ed_50 <- paste(all_genes_plus_all_chr_merge_50b$CNV_2ed_50bins,all_genes_plus_all_chr_merge_50b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_50b$CNV_3rd_50 <- paste(all_genes_plus_all_chr_merge_50b$CNV_3rd_50bins,all_genes_plus_all_chr_merge_50b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_50b$CNV_4th_50 <- paste(all_genes_plus_all_chr_merge_50b$CNV_4th_50bins,all_genes_plus_all_chr_merge_50b$chr_anno,sep="_")
  names <- c("CNV_1st_50","CNV_2ed_50","CNV_3rd_50","CNV_4th_50")
  all_genes_plus_all_chr_merge_50b_sel <- all_genes_plus_all_chr_merge_50b[,union(names,"hgnc_symbol")]

  rbind_all_chr_extract_counts <- function(name_num){
    sel_name <- names[name_num]
    CNV_segments_all <- unique(all_genes_plus_all_chr_merge_50b_sel[,sel_name])
    message(sel_name," is starting")
    sel_counts_all <- c()
    for (i in CNV_segments_all){
      sel_genes <- all_genes_plus_all_chr_merge_50b_sel[all_genes_plus_all_chr_merge_50b_sel[,name_num]==i,]$hgnc_symbol
        sel_genes <- intersect(sel_genes,rownames(epi_nor_data))
        if (length(sel_genes) != 1){
            sel_counts <- epi_nor_data[sel_genes,]
            sel_counts <- as.data.frame(colSums(sel_counts)/length(sel_genes))
            colnames(sel_counts) <- i
        }else{
          sel_counts <- as.data.frame(epi_nor_data[sel_genes,]/length(sel_genes))
          colnames(sel_counts) <- i
        }
        sel_counts <- as.data.frame(t(sel_counts))
        sel_counts_all <- rbind(sel_counts_all,sel_counts)
        #print(message(i," is done"))
    }
    message(sel_name," is done")
    return(sel_counts_all)
  }
  sel_counts_all <- future_lapply(1:length(names),rbind_all_chr_extract_counts)
  sel_counts_Seurat_obj <- do.call(rbind,sel_counts_all)
  sel_counts_Seurat_obj_50_bins <- sel_counts_Seurat_obj

  
  ##100 genes a bin
  all_chr_100b <- function(i){
    tmp <- subset(sel_gene_pos,chr_anno==i)
    tmp$order <- 1:nrow(tmp)
    rownames(tmp) <- tmp$ensembl_gene_id
    seed <- c(0,30,60,90)
    seed_names <- c("1st","2ed","3rd","4th")
    all_genes_plus_add_all <- c()
    for (order in 1:length(seed)){
      sel_order <- seed[order]
      if (nrow(tmp) > sel_order){
        if (sel_order != 0){
          sel_tmp <- tmp[1:sel_order,]
          sel_tmp$new_CNV <- paste("CNV",seed_names[order],"100bins_gap",sep="_")
          others_tmp <- tmp[setdiff(rownames(tmp),rownames(sel_tmp)),]
          }else{
            others_tmp <- tmp
            sel_tmp <- c()
        }
        }else{
          others_tmp <- tmp
          sel_tmp <- c()
      }
      if (nrow(others_tmp) > 100){
        tmp_order <- c()
        all_genes <- c()
        for(num in 1:floor(nrow(others_tmp)/100)){
          num_genes <- num*100
          order_genes <- rownames(others_tmp)[1:num_genes]
          order_genes <- setdiff(order_genes,tmp_order)
          tmp_genes <- others_tmp[order_genes,]
          tmp_genes$new_CNV <- paste("CNV",seed_names[order],"100bins",num,sep="_")
          tmp_order <- union(tmp_order,rownames(tmp_genes))
          all_genes <- rbind(all_genes,tmp_genes)
        }
        rownames(others_tmp) <- others_tmp$ensembl_gene_id
        if (length(setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id))==0){
          tmp_remind <- c()
          }else{
          tmp_remind <- others_tmp[setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id),]
          tmp_remind$new_CNV <- paste("CNV",seed_names[order],"100bins",round(nrow(others_tmp)/100),sep="_")
          }
        all_genes_plus <- rbind(all_genes,tmp_remind)
        all_genes_plus <- rbind(all_genes_plus,sel_tmp)
        }else{
          all_genes_plus <- others_tmp
          all_genes_plus$new_CNV <- paste("CNV",seed_names[order],"100bins",round(nrow(tmp)/100),sep="_")
          all_genes_plus <- rbind(all_genes_plus,sel_tmp)
      }
      rownames(all_genes_plus) <- all_genes_plus$ensembl_gene_id
      add_new <- data.frame(all_genes_plus[,ncol(all_genes_plus)],row.names=rownames(all_genes_plus))
      colnames(add_new) <- paste("CNV",seed_names[order],"100bins",sep="_")
      all_genes_plus <- cbind(all_genes_plus,add_new)
      all_genes_plus_sel <- as.data.frame(all_genes_plus[,setdiff(colnames(all_genes_plus),rownames(all_genes_plus_add_all))])
      if (ncol(all_genes_plus_sel)==1){
        colnames(all_genes_plus_sel) <- paste("CNV",seed_names[order],"100bins",sep="_")
        }else{
          all_genes_plus_sel <- all_genes_plus_sel
      }
      all_genes_plus_sel <- as.data.frame(t(all_genes_plus_sel))
      colnames(all_genes_plus_sel) <- rownames(all_genes_plus)
      all_genes_plus_add_all <- rbind(all_genes_plus_add_all,all_genes_plus_sel)
      #message(i,seed_names[order], " is done")
    }
    all_genes_plus_all_chr <- as.data.frame(t(all_genes_plus_add_all))
    return(all_genes_plus_all_chr)
    message(i," is done")
}

  all_chr_merge_100b <- future_lapply(unique(sel_gene_pos$chr_anno),all_chr_100b)
  all_genes_plus_all_chr_merge_100b <- do.call(rbind,all_chr_merge_100b)
  all_genes_plus_all_chr_merge_100b$CNV_1st_100 <- paste(all_genes_plus_all_chr_merge_100b$CNV_1st_100bins,all_genes_plus_all_chr_merge_100b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_100b$CNV_2ed_100 <- paste(all_genes_plus_all_chr_merge_100b$CNV_2ed_100bins,all_genes_plus_all_chr_merge_100b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_100b$CNV_3rd_100 <- paste(all_genes_plus_all_chr_merge_100b$CNV_3rd_100bins,all_genes_plus_all_chr_merge_100b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_100b$CNV_4th_100 <- paste(all_genes_plus_all_chr_merge_100b$CNV_4th_100bins,all_genes_plus_all_chr_merge_100b$chr_anno,sep="_")
  names <- c("CNV_1st_100","CNV_2ed_100","CNV_3rd_100","CNV_4th_100")
  all_genes_plus_all_chr_merge_100b_sel <- all_genes_plus_all_chr_merge_100b[,union(names,"hgnc_symbol")]

  rbind_all_chr_extract_counts <- function(name_num){
    sel_name <- names[name_num]
    CNV_segments_all <- unique(all_genes_plus_all_chr_merge_100b_sel[,sel_name])
    message(sel_name," is starting")
    sel_counts_all <- c()
    for (i in CNV_segments_all){
      sel_genes <- all_genes_plus_all_chr_merge_100b_sel[all_genes_plus_all_chr_merge_100b_sel[,name_num]==i,]$hgnc_symbol
        sel_genes <- intersect(sel_genes,rownames(epi_nor_data))
        if (length(sel_genes) != 1){
            sel_counts <- epi_nor_data[sel_genes,]
            sel_counts <- as.data.frame(colSums(sel_counts)/length(sel_genes))
            colnames(sel_counts) <- i
        }else{
          sel_counts <- as.data.frame(epi_nor_data[sel_genes,]/length(sel_genes))
          colnames(sel_counts) <- i
        }
        sel_counts <- as.data.frame(t(sel_counts))
        sel_counts_all <- rbind(sel_counts_all,sel_counts)
        #print(message(i," is done"))
    }
    message(sel_name," is done")
    return(sel_counts_all)
  }
  sel_counts_all <- future_lapply(1:length(names),rbind_all_chr_extract_counts)
  sel_counts_Seurat_obj <- do.call(rbind,sel_counts_all)
  sel_counts_Seurat_obj_100_bins <- sel_counts_Seurat_obj


  ##150 genes a bin
  all_chr_150b <- function(i){
    tmp <- subset(sel_gene_pos,chr_anno==i)
    tmp$order <- 1:nrow(tmp)
    rownames(tmp) <- tmp$ensembl_gene_id
    seed <- c(0,30,60,90)
    seed_names <- c("1st","2ed","3rd","4th")
    all_genes_plus_add_all <- c()
    for (order in 1:length(seed)){
      sel_order <- seed[order]
      if (nrow(tmp) > sel_order){
        if (sel_order != 0){
          sel_tmp <- tmp[1:sel_order,]
          sel_tmp$new_CNV <- paste("CNV",seed_names[order],"150bins_gap",sep="_")
          others_tmp <- tmp[setdiff(rownames(tmp),rownames(sel_tmp)),]
          }else{
            others_tmp <- tmp
            sel_tmp <- c()
        }
        }else{
          others_tmp <- tmp
          sel_tmp <- c()
      }
      if (nrow(others_tmp) > 150){
        tmp_order <- c()
        all_genes <- c()
        for(num in 1:floor(nrow(others_tmp)/150)){
          num_genes <- num*150
          order_genes <- rownames(others_tmp)[1:num_genes]
          order_genes <- setdiff(order_genes,tmp_order)
          tmp_genes <- others_tmp[order_genes,]
          tmp_genes$new_CNV <- paste("CNV",seed_names[order],"150bins",num,sep="_")
          tmp_order <- union(tmp_order,rownames(tmp_genes))
          all_genes <- rbind(all_genes,tmp_genes)
        }
        rownames(others_tmp) <- others_tmp$ensembl_gene_id
        if (length(setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id))==0){
          tmp_remind <- c()
          }else{
          tmp_remind <- others_tmp[setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id),]
          tmp_remind$new_CNV <- paste("CNV",seed_names[order],"150bins",round(nrow(others_tmp)/150),sep="_")
          }
        all_genes_plus <- rbind(all_genes,tmp_remind)
        all_genes_plus <- rbind(all_genes_plus,sel_tmp)
        }else{
          all_genes_plus <- others_tmp
          all_genes_plus$new_CNV <- paste("CNV",seed_names[order],"150bins",round(nrow(tmp)/150),sep="_")
          all_genes_plus <- rbind(all_genes_plus,sel_tmp)
      }
      rownames(all_genes_plus) <- all_genes_plus$ensembl_gene_id
      add_new <- data.frame(all_genes_plus[,ncol(all_genes_plus)],row.names=rownames(all_genes_plus))
      colnames(add_new) <- paste("CNV",seed_names[order],"150bins",sep="_")
      all_genes_plus <- cbind(all_genes_plus,add_new)
      all_genes_plus_sel <- as.data.frame(all_genes_plus[,setdiff(colnames(all_genes_plus),rownames(all_genes_plus_add_all))])
      if (ncol(all_genes_plus_sel)==1){
        colnames(all_genes_plus_sel) <- paste("CNV",seed_names[order],"150bins",sep="_")
        }else{
          all_genes_plus_sel <- all_genes_plus_sel
      }
      all_genes_plus_sel <- as.data.frame(t(all_genes_plus_sel))
      colnames(all_genes_plus_sel) <- rownames(all_genes_plus)
      all_genes_plus_add_all <- rbind(all_genes_plus_add_all,all_genes_plus_sel)
      #message(i,seed_names[order], " is done")
    }
    all_genes_plus_all_chr <- as.data.frame(t(all_genes_plus_add_all))
    return(all_genes_plus_all_chr)
    message(i," is done")
}

  all_chr_merge_150b <- future_lapply(unique(sel_gene_pos$chr_anno),all_chr_150b)
  all_genes_plus_all_chr_merge_150b <- do.call(rbind,all_chr_merge_150b)
  all_genes_plus_all_chr_merge_150b$CNV_1st_150 <- paste(all_genes_plus_all_chr_merge_150b$CNV_1st_150bins,all_genes_plus_all_chr_merge_150b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_150b$CNV_2ed_150 <- paste(all_genes_plus_all_chr_merge_150b$CNV_2ed_150bins,all_genes_plus_all_chr_merge_150b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_150b$CNV_3rd_150 <- paste(all_genes_plus_all_chr_merge_150b$CNV_3rd_150bins,all_genes_plus_all_chr_merge_150b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_150b$CNV_4th_150 <- paste(all_genes_plus_all_chr_merge_150b$CNV_4th_150bins,all_genes_plus_all_chr_merge_150b$chr_anno,sep="_")
  names <- c("CNV_1st_150","CNV_2ed_150","CNV_3rd_150","CNV_4th_150")
  all_genes_plus_all_chr_merge_150b_sel <- all_genes_plus_all_chr_merge_150b[,union(names,"hgnc_symbol")]

  rbind_all_chr_extract_counts <- function(name_num){
    sel_name <- names[name_num]
    CNV_segments_all <- unique(all_genes_plus_all_chr_merge_150b_sel[,sel_name])
    message(sel_name," is starting")
    sel_counts_all <- c()
    for (i in CNV_segments_all){
      sel_genes <- all_genes_plus_all_chr_merge_150b_sel[all_genes_plus_all_chr_merge_150b_sel[,name_num]==i,]$hgnc_symbol
        sel_genes <- intersect(sel_genes,rownames(epi_nor_data))
        if (length(sel_genes) != 1){
            sel_counts <- epi_nor_data[sel_genes,]
            sel_counts <- as.data.frame(colSums(sel_counts)/length(sel_genes))
            colnames(sel_counts) <- i
        }else{
          sel_counts <- as.data.frame(epi_nor_data[sel_genes,]/length(sel_genes))
          colnames(sel_counts) <- i
        }
        sel_counts <- as.data.frame(t(sel_counts))
        sel_counts_all <- rbind(sel_counts_all,sel_counts)
        #print(message(i," is done"))
    }
    message(sel_name," is done")
    return(sel_counts_all)
  }
  sel_counts_all <- future_lapply(1:length(names),rbind_all_chr_extract_counts)
  sel_counts_Seurat_obj <- do.call(rbind,sel_counts_all)
  sel_counts_Seurat_obj_150_bins <- sel_counts_Seurat_obj

  ##200 genes a bin
  all_chr_200b <- function(i){
    tmp <- subset(sel_gene_pos,chr_anno==i)
    tmp$order <- 1:nrow(tmp)
    rownames(tmp) <- tmp$ensembl_gene_id
    seed <- c(0,30,60,90)
    seed_names <- c("1st","2ed","3rd","4th")
    all_genes_plus_add_all <- c()
    for (order in 1:length(seed)){
      sel_order <- seed[order]
      if (nrow(tmp) > sel_order){
        if (sel_order != 0){
          sel_tmp <- tmp[1:sel_order,]
          sel_tmp$new_CNV <- paste("CNV",seed_names[order],"200bins_gap",sep="_")
          others_tmp <- tmp[setdiff(rownames(tmp),rownames(sel_tmp)),]
          }else{
            others_tmp <- tmp
            sel_tmp <- c()
        }
        }else{
          others_tmp <- tmp
          sel_tmp <- c()
      }
      if (nrow(others_tmp) > 200){
        tmp_order <- c()
        all_genes <- c()
        for(num in 1:floor(nrow(others_tmp)/200)){
          num_genes <- num*200
          order_genes <- rownames(others_tmp)[1:num_genes]
          order_genes <- setdiff(order_genes,tmp_order)
          tmp_genes <- others_tmp[order_genes,]
          tmp_genes$new_CNV <- paste("CNV",seed_names[order],"200bins",num,sep="_")
          tmp_order <- union(tmp_order,rownames(tmp_genes))
          all_genes <- rbind(all_genes,tmp_genes)
        }
        rownames(others_tmp) <- others_tmp$ensembl_gene_id
        if (length(setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id))==0){
          tmp_remind <- c()
          }else{
          tmp_remind <- others_tmp[setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id),]
          tmp_remind$new_CNV <- paste("CNV",seed_names[order],"200bins",round(nrow(others_tmp)/200),sep="_")
          }
        all_genes_plus <- rbind(all_genes,tmp_remind)
        all_genes_plus <- rbind(all_genes_plus,sel_tmp)
        }else{
          all_genes_plus <- others_tmp
          all_genes_plus$new_CNV <- paste("CNV",seed_names[order],"200bins",round(nrow(tmp)/200),sep="_")
          all_genes_plus <- rbind(all_genes_plus,sel_tmp)
      }
      rownames(all_genes_plus) <- all_genes_plus$ensembl_gene_id
      add_new <- data.frame(all_genes_plus[,ncol(all_genes_plus)],row.names=rownames(all_genes_plus))
      colnames(add_new) <- paste("CNV",seed_names[order],"200bins",sep="_")
      all_genes_plus <- cbind(all_genes_plus,add_new)
      all_genes_plus_sel <- as.data.frame(all_genes_plus[,setdiff(colnames(all_genes_plus),rownames(all_genes_plus_add_all))])
      if (ncol(all_genes_plus_sel)==1){
        colnames(all_genes_plus_sel) <- paste("CNV",seed_names[order],"200bins",sep="_")
        }else{
          all_genes_plus_sel <- all_genes_plus_sel
      }
      all_genes_plus_sel <- as.data.frame(t(all_genes_plus_sel))
      colnames(all_genes_plus_sel) <- rownames(all_genes_plus)
      all_genes_plus_add_all <- rbind(all_genes_plus_add_all,all_genes_plus_sel)
      #message(i,seed_names[order], " is done")
    }
    all_genes_plus_all_chr <- as.data.frame(t(all_genes_plus_add_all))
    return(all_genes_plus_all_chr)
    message(i," is done")
}

  all_chr_merge_200b <- future_lapply(unique(sel_gene_pos$chr_anno),all_chr_200b)
  all_genes_plus_all_chr_merge_200b <- do.call(rbind,all_chr_merge_200b)
  all_genes_plus_all_chr_merge_200b$CNV_1st_200 <- paste(all_genes_plus_all_chr_merge_200b$CNV_1st_200bins,all_genes_plus_all_chr_merge_200b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_200b$CNV_2ed_200 <- paste(all_genes_plus_all_chr_merge_200b$CNV_2ed_200bins,all_genes_plus_all_chr_merge_200b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_200b$CNV_3rd_200 <- paste(all_genes_plus_all_chr_merge_200b$CNV_3rd_200bins,all_genes_plus_all_chr_merge_200b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_200b$CNV_4th_200 <- paste(all_genes_plus_all_chr_merge_200b$CNV_4th_200bins,all_genes_plus_all_chr_merge_200b$chr_anno,sep="_")
  names <- c("CNV_1st_200","CNV_2ed_200","CNV_3rd_200","CNV_4th_200")
  all_genes_plus_all_chr_merge_200b_sel <- all_genes_plus_all_chr_merge_200b[,union(names,"hgnc_symbol")]

  rbind_all_chr_extract_counts <- function(name_num){
    sel_name <- names[name_num]
    CNV_segments_all <- unique(all_genes_plus_all_chr_merge_200b_sel[,sel_name])
    message(sel_name," is starting")
    sel_counts_all <- c()
    for (i in CNV_segments_all){
      sel_genes <- all_genes_plus_all_chr_merge_200b_sel[all_genes_plus_all_chr_merge_200b_sel[,name_num]==i,]$hgnc_symbol
        sel_genes <- intersect(sel_genes,rownames(epi_nor_data))
        if (length(sel_genes) != 1){
            sel_counts <- epi_nor_data[sel_genes,]
            sel_counts <- as.data.frame(colSums(sel_counts)/length(sel_genes))
            colnames(sel_counts) <- i
        }else{
          sel_counts <- as.data.frame(epi_nor_data[sel_genes,]/length(sel_genes))
          colnames(sel_counts) <- i
        }
        sel_counts <- as.data.frame(t(sel_counts))
        sel_counts_all <- rbind(sel_counts_all,sel_counts)
        #print(message(i," is done"))
    }
    message(sel_name," is done")
    return(sel_counts_all)
  }
  sel_counts_all <- future_lapply(1:length(names),rbind_all_chr_extract_counts)
  sel_counts_Seurat_obj <- do.call(rbind,sel_counts_all)
  sel_counts_Seurat_obj_200_bins <- sel_counts_Seurat_obj
  ##250 genes a bin
  all_chr_250b <- function(i){
    tmp <- subset(sel_gene_pos,chr_anno==i)
    tmp$order <- 1:nrow(tmp)
    rownames(tmp) <- tmp$ensembl_gene_id
    seed <- c(0,30,60,90)
    seed_names <- c("1st","2ed","3rd","4th")
    all_genes_plus_add_all <- c()
    for (order in 1:length(seed)){
      sel_order <- seed[order]
      if (nrow(tmp) > sel_order){
        if (sel_order != 0){
          sel_tmp <- tmp[1:sel_order,]
          sel_tmp$new_CNV <- paste("CNV",seed_names[order],"250bins_gap",sep="_")
          others_tmp <- tmp[setdiff(rownames(tmp),rownames(sel_tmp)),]
          }else{
            others_tmp <- tmp
            sel_tmp <- c()
        }
        }else{
          others_tmp <- tmp
          sel_tmp <- c()
      }
      if (nrow(others_tmp) > 250){
        tmp_order <- c()
        all_genes <- c()
        for(num in 1:floor(nrow(others_tmp)/250)){
          num_genes <- num*250
          order_genes <- rownames(others_tmp)[1:num_genes]
          order_genes <- setdiff(order_genes,tmp_order)
          tmp_genes <- others_tmp[order_genes,]
          tmp_genes$new_CNV <- paste("CNV",seed_names[order],"250bins",num,sep="_")
          tmp_order <- union(tmp_order,rownames(tmp_genes))
          all_genes <- rbind(all_genes,tmp_genes)
        }
        rownames(others_tmp) <- others_tmp$ensembl_gene_id
        if (length(setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id))==0){
          tmp_remind <- c()
          }else{
          tmp_remind <- others_tmp[setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id),]
          tmp_remind$new_CNV <- paste("CNV",seed_names[order],"250bins",round(nrow(others_tmp)/250),sep="_")
          }
        all_genes_plus <- rbind(all_genes,tmp_remind)
        all_genes_plus <- rbind(all_genes_plus,sel_tmp)
        }else{
          all_genes_plus <- others_tmp
          all_genes_plus$new_CNV <- paste("CNV",seed_names[order],"250bins",round(nrow(tmp)/250),sep="_")
          all_genes_plus <- rbind(all_genes_plus,sel_tmp)
      }
      rownames(all_genes_plus) <- all_genes_plus$ensembl_gene_id
      add_new <- data.frame(all_genes_plus[,ncol(all_genes_plus)],row.names=rownames(all_genes_plus))
      colnames(add_new) <- paste("CNV",seed_names[order],"250bins",sep="_")
      all_genes_plus <- cbind(all_genes_plus,add_new)
      all_genes_plus_sel <- as.data.frame(all_genes_plus[,setdiff(colnames(all_genes_plus),rownames(all_genes_plus_add_all))])
      if (ncol(all_genes_plus_sel)==1){
        colnames(all_genes_plus_sel) <- paste("CNV",seed_names[order],"250bins",sep="_")
        }else{
          all_genes_plus_sel <- all_genes_plus_sel
      }
      all_genes_plus_sel <- as.data.frame(t(all_genes_plus_sel))
      colnames(all_genes_plus_sel) <- rownames(all_genes_plus)
      all_genes_plus_add_all <- rbind(all_genes_plus_add_all,all_genes_plus_sel)
      #message(i,seed_names[order], " is done")
    }
    all_genes_plus_all_chr <- as.data.frame(t(all_genes_plus_add_all))
    return(all_genes_plus_all_chr)
    message(i," is done")
}

  all_chr_merge_250b <- future_lapply(unique(sel_gene_pos$chr_anno),all_chr_250b)
  all_genes_plus_all_chr_merge_250b <- do.call(rbind,all_chr_merge_250b)
  all_genes_plus_all_chr_merge_250b$CNV_1st_250 <- paste(all_genes_plus_all_chr_merge_250b$CNV_1st_250bins,all_genes_plus_all_chr_merge_250b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_250b$CNV_2ed_250 <- paste(all_genes_plus_all_chr_merge_250b$CNV_2ed_250bins,all_genes_plus_all_chr_merge_250b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_250b$CNV_3rd_250 <- paste(all_genes_plus_all_chr_merge_250b$CNV_3rd_250bins,all_genes_plus_all_chr_merge_250b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_250b$CNV_4th_250 <- paste(all_genes_plus_all_chr_merge_250b$CNV_4th_250bins,all_genes_plus_all_chr_merge_250b$chr_anno,sep="_")
  names <- c("CNV_1st_250","CNV_2ed_250","CNV_3rd_250","CNV_4th_250")
  all_genes_plus_all_chr_merge_250b_sel <- all_genes_plus_all_chr_merge_250b[,union(names,"hgnc_symbol")]

  rbind_all_chr_extract_counts <- function(name_num){
    sel_name <- names[name_num]
    CNV_segments_all <- unique(all_genes_plus_all_chr_merge_250b_sel[,sel_name])
    message(sel_name," is starting")
    sel_counts_all <- c()
    for (i in CNV_segments_all){
      sel_genes <- all_genes_plus_all_chr_merge_250b_sel[all_genes_plus_all_chr_merge_250b_sel[,name_num]==i,]$hgnc_symbol
        sel_genes <- intersect(sel_genes,rownames(epi_nor_data))
        if (length(sel_genes) != 1){
            sel_counts <- epi_nor_data[sel_genes,]
            sel_counts <- as.data.frame(colSums(sel_counts)/length(sel_genes))
            colnames(sel_counts) <- i
        }else{
          sel_counts <- as.data.frame(epi_nor_data[sel_genes,]/length(sel_genes))
          colnames(sel_counts) <- i
        }
        sel_counts <- as.data.frame(t(sel_counts))
        sel_counts_all <- rbind(sel_counts_all,sel_counts)
        #print(message(i," is done"))
    }
    message(sel_name," is done")
    return(sel_counts_all)
  }
  sel_counts_all <- future_lapply(1:length(names),rbind_all_chr_extract_counts)
  sel_counts_Seurat_obj <- do.call(rbind,sel_counts_all)
  sel_counts_Seurat_obj_250_bins <- sel_counts_Seurat_obj
  ##300 genes a bin
  all_chr_300b <- function(i){
    tmp <- subset(sel_gene_pos,chr_anno==i)
    tmp$order <- 1:nrow(tmp)
    rownames(tmp) <- tmp$ensembl_gene_id
    seed <- c(0,30,60,90)
    seed_names <- c("1st","2ed","3rd","4th")
    all_genes_plus_add_all <- c()
    for (order in 1:length(seed)){
      sel_order <- seed[order]
      if (nrow(tmp) > sel_order){
        if (sel_order != 0){
          sel_tmp <- tmp[1:sel_order,]
          sel_tmp$new_CNV <- paste("CNV",seed_names[order],"300bins_gap",sep="_")
          others_tmp <- tmp[setdiff(rownames(tmp),rownames(sel_tmp)),]
          }else{
            others_tmp <- tmp
            sel_tmp <- c()
        }
        }else{
          others_tmp <- tmp
          sel_tmp <- c()
      }
      if (nrow(others_tmp) > 300){
        tmp_order <- c()
        all_genes <- c()
        for(num in 1:floor(nrow(others_tmp)/300)){
          num_genes <- num*300
          order_genes <- rownames(others_tmp)[1:num_genes]
          order_genes <- setdiff(order_genes,tmp_order)
          tmp_genes <- others_tmp[order_genes,]
          tmp_genes$new_CNV <- paste("CNV",seed_names[order],"300bins",num,sep="_")
          tmp_order <- union(tmp_order,rownames(tmp_genes))
          all_genes <- rbind(all_genes,tmp_genes)
        }
        rownames(others_tmp) <- others_tmp$ensembl_gene_id
        if (length(setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id))==0){
          tmp_remind <- c()
          }else{
          tmp_remind <- others_tmp[setdiff(others_tmp$ensembl_gene_id,all_genes$ensembl_gene_id),]
          tmp_remind$new_CNV <- paste("CNV",seed_names[order],"300bins",round(nrow(others_tmp)/300),sep="_")
          }
        all_genes_plus <- rbind(all_genes,tmp_remind)
        all_genes_plus <- rbind(all_genes_plus,sel_tmp)
        }else{
          all_genes_plus <- others_tmp
          all_genes_plus$new_CNV <- paste("CNV",seed_names[order],"300bins",round(nrow(tmp)/300),sep="_")
          all_genes_plus <- rbind(all_genes_plus,sel_tmp)
      }
      rownames(all_genes_plus) <- all_genes_plus$ensembl_gene_id
      add_new <- data.frame(all_genes_plus[,ncol(all_genes_plus)],row.names=rownames(all_genes_plus))
      colnames(add_new) <- paste("CNV",seed_names[order],"300bins",sep="_")
      all_genes_plus <- cbind(all_genes_plus,add_new)
      all_genes_plus_sel <- as.data.frame(all_genes_plus[,setdiff(colnames(all_genes_plus),rownames(all_genes_plus_add_all))])
      if (ncol(all_genes_plus_sel)==1){
        colnames(all_genes_plus_sel) <- paste("CNV",seed_names[order],"300bins",sep="_")
        }else{
          all_genes_plus_sel <- all_genes_plus_sel
      }
      all_genes_plus_sel <- as.data.frame(t(all_genes_plus_sel))
      colnames(all_genes_plus_sel) <- rownames(all_genes_plus)
      all_genes_plus_add_all <- rbind(all_genes_plus_add_all,all_genes_plus_sel)
      #message(i,seed_names[order], " is done")
    }
    all_genes_plus_all_chr <- as.data.frame(t(all_genes_plus_add_all))
    return(all_genes_plus_all_chr)
    message(i," is done")
}

  all_chr_merge_300b <- future_lapply(unique(sel_gene_pos$chr_anno),all_chr_300b)
  all_genes_plus_all_chr_merge_300b <- do.call(rbind,all_chr_merge_300b)
  all_genes_plus_all_chr_merge_300b$CNV_1st_300 <- paste(all_genes_plus_all_chr_merge_300b$CNV_1st_300bins,all_genes_plus_all_chr_merge_300b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_300b$CNV_2ed_300 <- paste(all_genes_plus_all_chr_merge_300b$CNV_2ed_300bins,all_genes_plus_all_chr_merge_300b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_300b$CNV_3rd_300 <- paste(all_genes_plus_all_chr_merge_300b$CNV_3rd_300bins,all_genes_plus_all_chr_merge_300b$chr_anno,sep="_")
  all_genes_plus_all_chr_merge_300b$CNV_4th_300 <- paste(all_genes_plus_all_chr_merge_300b$CNV_4th_300bins,all_genes_plus_all_chr_merge_300b$chr_anno,sep="_")
  names <- c("CNV_1st_300","CNV_2ed_300","CNV_3rd_300","CNV_4th_300")
  all_genes_plus_all_chr_merge_300b_sel <- all_genes_plus_all_chr_merge_300b[,union(names,"hgnc_symbol")]

  rbind_all_chr_extract_counts <- function(name_num){
    sel_name <- names[name_num]
    CNV_segments_all <- unique(all_genes_plus_all_chr_merge_300b_sel[,sel_name])
    message(sel_name," is starting")
    sel_counts_all <- c()
    for (i in CNV_segments_all){
      sel_genes <- all_genes_plus_all_chr_merge_300b_sel[all_genes_plus_all_chr_merge_300b_sel[,name_num]==i,]$hgnc_symbol
        sel_genes <- intersect(sel_genes,rownames(epi_nor_data))
        if (length(sel_genes) != 1){
            sel_counts <- epi_nor_data[sel_genes,]
            sel_counts <- as.data.frame(colSums(sel_counts)/length(sel_genes))
            colnames(sel_counts) <- i
        }else{
          sel_counts <- as.data.frame(epi_nor_data[sel_genes,]/length(sel_genes))
          colnames(sel_counts) <- i
        }
        sel_counts <- as.data.frame(t(sel_counts))
        sel_counts_all <- rbind(sel_counts_all,sel_counts)
        #print(message(i," is done"))
    }
    message(sel_name," is done")
    return(sel_counts_all)
  }
  sel_counts_all <- future_lapply(1:length(names),rbind_all_chr_extract_counts)
  sel_counts_Seurat_obj <- do.call(rbind,sel_counts_all)
  sel_counts_Seurat_obj_300_bins <- sel_counts_Seurat_obj

  sel_counts_Seurat_obj <- list(sel_counts_Seurat_obj_50_bins,sel_counts_Seurat_obj_100_bins,sel_counts_Seurat_obj_150_bins,
    sel_counts_Seurat_obj_200_bins,sel_counts_Seurat_obj_250_bins,sel_counts_Seurat_obj_300_bins)
  sel_counts_Seurat_obj <- do.call(rbind,sel_counts_Seurat_obj)
  meta_info <- Seurat_obj@meta.data
  meta_info <- meta_info[colnames(sel_counts_Seurat_obj),]
  Pseudo_CNV_obj <- CreateSeuratObject(counts = sel_counts_Seurat_obj, meta.data=meta_info,project = project)
  if (anno_pos_need==TRUE){
    all_peuso_CNV <- list(all_genes_plus_all_chr_merge_100b,all_genes_plus_all_chr_merge_150b,all_genes_plus_all_chr_merge_200b,
      all_genes_plus_all_chr_merge_250b,all_genes_plus_all_chr_merge_300b)
    cbind_processed <- function(x){
      x <- x[,-c(1:14)]
      return(x)
    }
    all_peuso_CNV <- future_lapply(all_peuso_CNV,cbind_processed)
    all_genes_plus_all_chr_merge_50b <- all_genes_plus_all_chr_merge_50b[,-c(10:14)]
    all_peuso_CNV <- c(list(all_genes_plus_all_chr_merge_50b),all_peuso_CNV)
    all_peuso_CNV_list <- do.call(cbind,all_peuso_CNV)
    mcsaveRDS(all_peuso_CNV_list,paste0(save_path,"/",file_prefix,"peuso_CNV_list.rds"),mc.cores=20)
  }
  return(Pseudo_CNV_obj)
}







add_CNV_ASSAY <- function(RNA_Seurat_obj=RNA_Seurat_obj,CNV_Seurat_obj=CNV_Seurat_obj,RBE=RBE,dims=dims)
{
	require(Seurat)
	require(future.apply)
	if (RBE==TRUE){
		if (ncol(RNA_Seurat_obj)==ncol(CNV_Seurat_obj)){
			Seurat_CNV_RNA <- RNA_Seurat_obj
			Seurat_CNV_RNA[["CNV_RD"]] <- CreateAssayObject(counts=as.matrix(GetAssayData(CNV_Seurat_obj,slot="counts")))
			Seurat_CNV_RNA <- Seurat_CNV_RNA %>%
			Seurat::NormalizeData(verbose = FALSE, assay = "CNV_RD") %>%
			FindVariableFeatures(selection.method = "vst", nfeatures = length(rownames(Seurat_CNV_RNA)),assay="CNV_RD") %>% 
		    ScaleData(verbose = TRUE, vars.to.regress = "nCount_CNV_RD", assay = "CNV_RD") %>% 
		    RunPCA(npcs = 30, verbose = TRUE, assay = "CNV_RD") %>%
		    RunHarmony("group", plot_convergence = TRUE, assay = "CNV_RD") %>% 
		    RunUMAP(reduction = "harmony", dims = dims, assay = "CNV_RD",reduction.name="CNV_harmony_umap") %>% 
		    RunTSNE(reduction = "harmony", dims = dims, assay = "CNV_RD",reduction.name="CNV_harmony_tsne") 
			DefaultAssay(Seurat_CNV_RNA) <- "CNV_RD"
			} else {
				stop("Cannot match same cells in CNV obj and RNA obj")
			}
			return(Seurat_CNV_RNA)
		} else {
			if (ncol(RNA_Seurat_obj)==ncol(CNV_Seurat_obj)){
				Seurat_CNV_RNA <- RNA_Seurat_obj
				Seurat_CNV_RNA[["CNV_RD"]] <- CreateAssayObject(counts=as.matrix(GetAssayData(CNV_Seurat_obj,slot="counts")))
				Seurat_CNV_RNA <- Seurat_CNV_RNA %>%
				Seurat::NormalizeData(verbose = FALSE, assay = "CNV_RD") %>%
				FindVariableFeatures(selection.method = "vst", nfeatures = length(rownames(Seurat_CNV_RNA)),assay="CNV_RD") %>% 
			    ScaleData(verbose = TRUE, vars.to.regress = "nCount_CNV_RD", assay = "CNV_RD") %>% 
			    RunPCA(npcs = 30, verbose = TRUE, assay = "CNV_RD") %>% 
				RunUMAP(dims = dims,reduction.name = "CNV_RD_UMAP", reduction.key = "CNV_RD_UMAP_", assay = "CNV_RD")%>% 
				RunTSNE(dims = dims,reduction.name = "CNV_RD_TSNE", reduction.key = "CNV_RD_TSNE_", assay = "CNV_RD")
				DefaultAssay(Seurat_CNV_RNA) <- "CNV_RD"
				} else {
					stop("Cannot match same cells in CNV obj and RNA obj")
				}
				return(Seurat_CNV_RNA)
			}
		}





DimPlot_FOR_Mut <- function (object, dims = c(1, 2), cells = NULL, cols = NULL,
    pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL,
    shape.by = NULL, order = NULL, label = FALSE, label.size = 4,
    repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26",
    sizes.highlight = 1, na.value = "grey50", combine = TRUE,
    ncol = NULL, ...)
{
    CheckDots(..., fxns = "CombinePlots")
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    reduction <- reduction %||% DefaultDimReduc(object = object)
    cells <- cells %||% colnames(x = object)
    data <- Embeddings(object = object[[reduction]])[cells, dims]
    data <- as.data.frame(x = data)
    dims <- paste0(Key(object = object[[reduction]]), dims)
    object[["ident"]] <- Idents(object = object)
    group.by <- group.by %||% "ident"
    data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
    for (group in group.by) {
        if (!is.factor(x = data[, group])) {
            data[, group] <- factor(x = data[, group])
        }
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    if (!is.null(x = split.by)) {
        data[, split.by] <- object[[split.by, drop = TRUE]]
    }
    levels_order <- c("No Call","ref/ref","alt/ref","alt/alt")
    levels_order <- levels_order[levels_order%in%levels(data[,3])]
    data[,3] <- factor(data[,3],levels=levels_order)
    data <- data[order(data[,3],decreasing=F),]
    plots <- lapply(X = group.by, FUN = function(x) {
        plot <- SingleDimPlot(data = data[, c(dims, x, split.by,
            shape.by)], dims = dims, col.by = x, cols = cols,
            pt.size = pt.size, shape.by = shape.by, order = order,
            label = FALSE, cells.highlight = cells.highlight,
            cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
            na.value = na.value) 
        if (label) {
            plot <- LabelClusters(plot = plot, id = x, repel = repel,
                size = label.size, split.by = split.by)
        }
        if (!is.null(x = split.by)) {
            plot <- plot + FacetTheme() + facet_wrap(facets = vars(!!sym(x = split.by)),
                ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                  length(x = unique(x = data[, split.by]))
                }
                else {
                  ncol
                })
        }
        return(plot)
    })
    if (combine) {
        plots <- CombinePlots(plots = plots, ncol = if (!is.null(x = split.by) &&
            length(x = group.by) > 1) {
            1
        }
        else {
            ncol
        }, ...)
    }
    return(plots)
}
environment(DimPlot_FOR_Mut) <- asNamespace('Seurat')

extract_Mut_info <- function(snv_seurat=snv_seurat,rna_seurat=rna_seurat,sel_genes=sel_genes,gene_info=gene_info,transfer_label=transfer_label){
	require(Seurat)
	require(future.apply)
	require(stringr)
	if (length(sel_genes)==1){
		counts_sel <- data.frame(GetAssayData(snv_seurat,slot="counts")[sel_genes,])
		colnames(counts_sel) <- sel_genes
		counts_sel <- as.data.frame(t(counts_sel))
		rownames(counts_sel) <- paste0(gene_info[sel_genes,]$SYMBOL,"_",rownames(counts_sel))

		} else {
		counts_sel <- as.data.frame(GetAssayData(snv_seurat,slot="counts")[sel_genes,])
		rownames(counts_sel) <- paste0(gene_info[sel_genes,]$SYMBOL,"_",rownames(counts_sel))
		}
	
	if (transfer_label==TRUE){
		counts_sel_trans <- future_apply(counts_sel,1,function(x){
		new_name <- str_replace(as.character(x), "0", "No Call")
		new_name <- str_replace(as.character(new_name), "1", "ref/ref")
		new_name <- str_replace(as.character(new_name), "2", "alt/alt")
		new_name <- str_replace(as.character(new_name), "3", "alt/ref")
		return(new_name)
		})
	} else {
		counts_sel_trans <- as.data.frame(t(counts_sel))
	}
	rownames(counts_sel_trans) <- colnames(counts_sel)
	if (length(sel_genes)==1){
		tmp <- as.data.frame(counts_sel_trans[rownames(rna_seurat@meta.data),])
		colnames(tmp) <- colnames(counts_sel_trans)
		rna_seurat@meta.data <- data.frame(tmp,row.names=rownames(rna_seurat@meta.data))
		} else {
			counts_sel_trans <- counts_sel_trans[rownames(rna_seurat@meta.data),]
			rna_seurat@meta.data <- data.frame(counts_sel_trans,row.names=rownames(rna_seurat@meta.data))
		}
	return(rna_seurat)
}





Visulize_Mut_Cluster <- function(snv_seurat=snv_seurat,rna_seurat=rna_seurat,reduction=reduction,
	sel_genes=sel_genes,gene_info=gene_info,sel_clu=sel_clu,split_T=split_T){
	require(pheatmap)
	require(ComplexHeatmap)
	require(BuenColors)
	require(scales)
	seurat_tmp <- rna_seurat
	extract_Clu_Mut <- function(i){
		sel_clu_num <- as.character(unique(gene_info[,sel_clu]))[i]
		sel_data <- gene_info[gene_info[,sel_clu]==sel_clu_num,]
		tmp_obj <- extract_Mut_info(snv_seurat=snv_seurat,
			rna_seurat=rna_seurat,
			sel_genes=sel_data$gene_pos,
			gene_info=sel_data,
			transfer_label=FALSE)
		meta_info <- tmp_obj@meta.data
		if (ncol(meta_info)==1){
			meta_info_sub <- as.data.frame(meta_info)
			colnames(meta_info_sub) <- colnames(tmp_obj@meta.data) 
			} else {
				meta_info_sub <- future_apply(meta_info,1,mean)
				meta_info_sub <- as.data.frame(meta_info_sub)
				meta_info_sub[meta_info_sub>0] <- 2
				meta_info_sub[meta_info_sub<1] <- 0
				colnames(meta_info_sub) <- sel_clu_num
				meta_info_sub$names <- rownames(meta_info_sub)
				tmp_obj@meta.data[[sel_clu_num]] <- meta_info_sub[rownames(tmp_obj@meta.data),1]
				meta_info_sub <- tmp_obj@meta.data
			}
			return(meta_info_sub)
	}
	meta_all <- future_lapply(1:length(unique(gene_info[,sel_clu])),extract_Clu_Mut)
	names(meta_all) <- as.character(unique(gene_info[,sel_clu]))
	meta_all_trans <- future_lapply(1:length(meta_all),function(x){
		tmp <- future_apply(meta_all[[x]],2,function(y){
			new_name <- str_replace(as.character(y), "0", "No Call")
			new_name <- str_replace(as.character(new_name), "2", "alt/alt")
			return(new_name)})
		rownames(tmp) <- rownames(meta_all[[x]])
		return(tmp)
		})
	names(meta_all_trans) <- as.character(unique(gene_info[,sel_clu]))
	col_sel <- c("lightgrey",hue_pal()(3)[3:1])
	names(col_sel) <- c("No Call","ref/ref","alt/ref","alt/alt")
	message("data processed done")
	message("printing begin")
	plots_all <- list()
		for (i in 1:length(meta_all_trans)){
			meta_all_trans_Sel <- meta_all_trans[[i]]
			if (ncol(meta_all_trans_Sel)==1){
				seurat_tmp@meta.data <- data.frame(meta_all_trans_Sel[rownames(seurat_tmp@meta.data),],row.names=rownames(seurat_tmp@meta.data))
				seurat_tmp@meta.data <- data.frame(seurat_tmp@meta.data,seurat_tmp@meta.data)
				colnames(seurat_tmp@meta.data) <- c(colnames(meta_all_trans_Sel),names(meta_all_trans)[i])
				colnames(seurat_tmp@meta.data) <- gsub(":",".",colnames(seurat_tmp@meta.data))
				} else {
					seurat_tmp@meta.data <- as.data.frame(meta_all_trans_Sel[rownames(seurat_tmp@meta.data),])
				}
			if(split_T==TRUE){
			plots <- list()
			for (j in 1:ncol(seurat_tmp@meta.data)){
				plots[[j]] <- DimPlot_FOR_Mut(object = seurat_tmp, group.by=colnames(seurat_tmp@meta.data)[j],reduction = reduction,pt.size=1,label = TRUE,repel=T,cols=col_sel) + labs(title=colnames(seurat_tmp@meta.data)[j])
				#message(names(meta_all_trans)[i]," ",colnames(seurat_tmp@meta.data)[j], " is done")
			}
			} else {
				pp <- DimPlot_FOR_Mut(object = seurat_tmp, group.by=names(meta_all_trans)[i],reduction = reduction,pt.size=1,label = TRUE,repel=T,cols=col_sel) + labs(title=names(meta_all_trans)[i])
				plots <- pp
			}
		plots_all[[i]] <- plots
		names(plots_all)[i] <- names(meta_all_trans)[i]
		message(names(meta_all_trans)[i], " is done")
	}
	return(plots_all)
}


