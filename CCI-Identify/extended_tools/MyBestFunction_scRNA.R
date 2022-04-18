XY_FeaturePlot <- function (object, features, dims = c(1, 2), cells = NULL, cols = c("lightgrey",
      "blue"), pt.size = NULL, min.cutoff = NA, max.cutoff = NA,
      reduction = NULL, split.by = NULL, shape.by = NULL, blend = FALSE,
      blend.threshold = 0.5, order = NULL, label = FALSE, label.size = 4,
      ncol = NULL, combine = TRUE, coord.fixed = FALSE, ...)
 {
      no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",
              size = 14, margin = margin(r = 7)))
      if (is.null(reduction)) {
          default_order <- c("umap", "tsne", "pca","dm")
          reducs <- which(default_order %in% names(object@reductions))
          reduction <- default_order[reducs[1]]
      }
      if (length(x = dims) != 2 || !is.numeric(x = dims)) {
          stop("'dims' must be a two-length integer vector")
      }
      if (blend && length(x = features) != 2) {
          stop("Blending feature plots only works with two features")
      }
      dims <- paste0(Key(object = object[[reduction]]), dims)
      cells <- cells %||% colnames(x = object)
      data <- FetchData(object = object, vars = c(dims, features),
          cells = cells)
      features <- colnames(x = data)[3:ncol(x = data)]
      min.cutoff <- mapply(FUN = function(cutoff, feature) {
          return(ifelse(test = is.na(x = cutoff), yes = min(data[,
              feature]), no = cutoff))
      }, cutoff = min.cutoff, feature = features)
      max.cutoff <- mapply(FUN = function(cutoff, feature) {
          return(ifelse(test = is.na(x = cutoff), yes = max(data[,
              feature]), no = cutoff))
      }, cutoff = max.cutoff, feature = features)
      check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
          max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
      if (length(x = check.lengths) != 1) {
          stop("There must be the same number of minimum and maximum cuttoffs as there are features")
      }
      brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
          ]$maxcolors, no = length(x = cols))
      data[, 3:ncol(x = data)] <- sapply(X = 3:ncol(x = data),
          FUN = function(index) {
              data.feature <- as.vector(x = data[, index])
              min.use <- SetQuantile(cutoff = min.cutoff[index -
                  2], data.feature)
              max.use <- SetQuantile(cutoff = max.cutoff[index -
                  2], data.feature)
              data.feature[data.feature < min.use] <- min.use
              data.feature[data.feature > max.use] <- max.use
              if (brewer.gran == 2) {
                  return(data.feature)
              }
              data.cut <- if (all(data.feature == 0)) {
                  0
              }
              else {
                  as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                    breaks = brewer.gran)))
              }
              return(data.cut)
          })
      colnames(x = data)[3:ncol(x = data)] <- features
      rownames(x = data) <- cells
      data$split <- if (is.null(x = split.by)) {
          RandomName()
      }
      else {
          switch(EXPR = split.by, ident = Idents(object = object)[cells],
              object[[split.by, drop = TRUE]][cells])
      }
      if (!is.factor(x = data$split)) {
          data$split <- factor(x = data$split)
      }
      plots <- vector(mode = "list", length = ifelse(test = blend,
          yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
      xlims <- c(min(data[, dims[1]]), max(data[,dims[1]]))
      ylims <- c(min(data[, dims[2]]), max(data[,dims[2]]))
      if (blend) {
          ncol <- 4
          color.matrix <- BlendMatrix(col.threshold = blend.threshold)
          colors <- list(color.matrix[, 1], color.matrix[1, ],
              as.vector(x = color.matrix))
      }
      for (i in 1:length(x = levels(x = data$split))) {
          ident <- levels(x = data$split)[i]
          data.plot <- data[as.character(x = data$split) == ident,
              , drop = FALSE]
          if (blend) {
              data.plot <- cbind(data.plot[, dims], BlendExpression(data = data.plot[,
                  features[1:2]]))
              features <- colnames(x = data.plot)[3:ncol(x = data.plot)]
          }
          for (j in 1:length(x = features)) {
              feature <- features[j]
              if (blend) {
                  cols.use <- as.numeric(x = as.character(x = data.plot[,
                    feature])) + 1
                  cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
              }
              else {
                  cols.use <- NULL
              }
              data.plot <- data.plot[order(data.plot[,j+2],decreasing=F),]
              plot <- SingleDimPlot(data = data.plot[, c(dims,
                  feature)], dims = dims, col.by = feature, pt.size = pt.size,
                  cols = cols.use, label.size = label.size) +
                  scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
                  theme_cowplot()
              if (length(x = levels(x = data$split)) > 1) {
                  plot <- plot + theme(panel.border = element_rect(fill = NA,
                    colour = "black"))
                  plot <- plot + if (i == 1) {
                    labs(title = feature)
                  }
                  else {
                    labs(title = NULL)
                  }
                  if (j == length(x = features) && !blend) {
                    suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident)) +
                      no.right)
                  }
                  if (j != 1) {
                    plot <- plot + theme(axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                      axis.title.y.left = element_blank())
                  }
                  if (i != length(x = levels(x = data$split))) {
                    plot <- plot + theme(axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                      axis.title.x = element_blank())
                  }
              }
              else {
                  plot <- plot + labs(title = feature)
              }
              if (!blend) {
                  plot <- plot + guides(color = NULL)
                  if (length(x = cols) == 1) {
                    plot <- plot + scale_color_brewer(palette = cols)
                  }
                  else if (length(x = cols) > 1) {
                    plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols,
                      guide = "colorbar"))
                  }
              }
              if (coord.fixed) {
                  plot <- plot + coord_fixed()
              }
              plot <- plot
              plots[[(length(x = features) * (i - 1)) + j]] <- plot
          }
      }
      if (blend) {
          blend.legend <- BlendMap(color.matrix = color.matrix)
          for (i in 1:length(x = levels(x = data$split))) {
              suppressMessages(expr = plots <- append(x = plots,
                  values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                    1, yes = levels(x = data$split)[i], no = "")),
                    expand = c(0, 0)) + labs(x = features[1], y = features[2],
                    title = if (i == 1) {
                      paste("Color threshold:", blend.threshold)
                    } else {
                      NULL
                    }) + no.right), after = 4 * i - 1))
          }
      }
      plots <- Filter(f = Negate(f = is.null), x = plots)
      if (combine) {
          if (is.null(x = ncol)) {
              ncol <- 2
              if (length(x = features) == 1) {
                  ncol <- 1
              }
              if (length(x = features) > 6) {
                  ncol <- 3
              }
              if (length(x = features) > 9) {
                  ncol <- 4
              }
          }
          ncol <- ifelse(test = is.null(x = split.by) || blend,
              yes = ncol, no = length(x = features))
          legend <- if (blend) {
              "none"
          }
          else {
              split.by %iff% "none"
          }
          plots <- CombinePlots(plots = plots, ncol = ncol, legend = legend,
              nrow = split.by %iff% length(x = levels(x = data$split)))
      }
      return(plots)
}
environment(XY_FeaturePlot) <- asNamespace('Seurat')

XY_FeaturePlot_no_order <- function (object, features, dims = c(1, 2), cells = NULL, cols = c("lightgrey",
      "blue"), pt.size = NULL, min.cutoff = NA, max.cutoff = NA,
      reduction = NULL, split.by = NULL, shape.by = NULL, blend = FALSE,
      blend.threshold = 0.5, order = NULL, label = FALSE, label.size = 4,
      ncol = NULL, combine = TRUE, coord.fixed = FALSE, ...)
 {
      no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",
              size = 14, margin = margin(r = 7)))
      if (is.null(reduction)) {
          default_order <- c("umap", "tsne", "pca","dm")
          reducs <- which(default_order %in% names(object@reductions))
          reduction <- default_order[reducs[1]]
      }
      if (length(x = dims) != 2 || !is.numeric(x = dims)) {
          stop("'dims' must be a two-length integer vector")
      }
      if (blend && length(x = features) != 2) {
          stop("Blending feature plots only works with two features")
      }
      dims <- paste0(Key(object = object[[reduction]]), dims)
      cells <- cells %||% colnames(x = object)
      data <- FetchData(object = object, vars = c(dims, features),
          cells = cells)
      features <- colnames(x = data)[3:ncol(x = data)]
      min.cutoff <- mapply(FUN = function(cutoff, feature) {
          return(ifelse(test = is.na(x = cutoff), yes = min(data[,
              feature]), no = cutoff))
      }, cutoff = min.cutoff, feature = features)
      max.cutoff <- mapply(FUN = function(cutoff, feature) {
          return(ifelse(test = is.na(x = cutoff), yes = max(data[,
              feature]), no = cutoff))
      }, cutoff = max.cutoff, feature = features)
      check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
          max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
      if (length(x = check.lengths) != 1) {
          stop("There must be the same number of minimum and maximum cuttoffs as there are features")
      }
      brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
          ]$maxcolors, no = length(x = cols))
      data[, 3:ncol(x = data)] <- sapply(X = 3:ncol(x = data),
          FUN = function(index) {
              data.feature <- as.vector(x = data[, index])
              min.use <- SetQuantile(cutoff = min.cutoff[index -
                  2], data.feature)
              max.use <- SetQuantile(cutoff = max.cutoff[index -
                  2], data.feature)
              data.feature[data.feature < min.use] <- min.use
              data.feature[data.feature > max.use] <- max.use
              if (brewer.gran == 2) {
                  return(data.feature)
              }
              data.cut <- if (all(data.feature == 0)) {
                  0
              }
              else {
                  as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                    breaks = brewer.gran)))
              }
              return(data.cut)
          })
      colnames(x = data)[3:ncol(x = data)] <- features
      rownames(x = data) <- cells
      data$split <- if (is.null(x = split.by)) {
          RandomName()
      }
      else {
          switch(EXPR = split.by, ident = Idents(object = object)[cells],
              object[[split.by, drop = TRUE]][cells])
      }
      if (!is.factor(x = data$split)) {
          data$split <- factor(x = data$split)
      }
      plots <- vector(mode = "list", length = ifelse(test = blend,
          yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
      xlims <- c(min(data[, dims[1]]), max(data[,dims[1]]))
      ylims <- c(min(data[, dims[2]]), max(data[,dims[2]]))
      if (blend) {
          ncol <- 4
          color.matrix <- BlendMatrix(col.threshold = blend.threshold)
          colors <- list(color.matrix[, 1], color.matrix[1, ],
              as.vector(x = color.matrix))
      }
      for (i in 1:length(x = levels(x = data$split))) {
          ident <- levels(x = data$split)[i]
          data.plot <- data[as.character(x = data$split) == ident,
              , drop = FALSE]
          if (blend) {
              data.plot <- cbind(data.plot[, dims], BlendExpression(data = data.plot[,
                  features[1:2]]))
              features <- colnames(x = data.plot)[3:ncol(x = data.plot)]
          }
          for (j in 1:length(x = features)) {
              feature <- features[j]
              if (blend) {
                  cols.use <- as.numeric(x = as.character(x = data.plot[,
                    feature])) + 1
                  cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
              }
              else {
                  cols.use <- NULL
              }
#              data.plot <- data.plot[order(data.plot[,j+2],decreasing=F),]
              plot <- SingleDimPlot(data = data.plot[, c(dims,
                  feature)], dims = dims, col.by = feature, pt.size = pt.size,
                  cols = cols.use, label.size = label.size) +
                  scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
                  theme_cowplot()
              if (length(x = levels(x = data$split)) > 1) {
                  plot <- plot + theme(panel.border = element_rect(fill = NA,
                    colour = "black"))
                  plot <- plot + if (i == 1) {
                    labs(title = feature)
                  }
                  else {
                    labs(title = NULL)
                  }
                  if (j == length(x = features) && !blend) {
                    suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident)) +
                      no.right)
                  }
                  if (j != 1) {
                    plot <- plot + theme(axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                      axis.title.y.left = element_blank())
                  }
                  if (i != length(x = levels(x = data$split))) {
                    plot <- plot + theme(axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                      axis.title.x = element_blank())
                  }
              }
              else {
                  plot <- plot + labs(title = feature)
              }
              if (!blend) {
                  plot <- plot + guides(color = NULL)
                  if (length(x = cols) == 1) {
                    plot <- plot + scale_color_brewer(palette = cols)
                  }
                  else if (length(x = cols) > 1) {
                    plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols,
                      guide = "colorbar"))
                  }
              }
              if (coord.fixed) {
                  plot <- plot + coord_fixed()
              }
              plot <- plot
              plots[[(length(x = features) * (i - 1)) + j]] <- plot
          }
      }
      if (blend) {
          blend.legend <- BlendMap(color.matrix = color.matrix)
          for (i in 1:length(x = levels(x = data$split))) {
              suppressMessages(expr = plots <- append(x = plots,
                  values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                    1, yes = levels(x = data$split)[i], no = "")),
                    expand = c(0, 0)) + labs(x = features[1], y = features[2],
                    title = if (i == 1) {
                      paste("Color threshold:", blend.threshold)
                    } else {
                      NULL
                    }) + no.right), after = 4 * i - 1))
          }
      }
      plots <- Filter(f = Negate(f = is.null), x = plots)
      if (combine) {
          if (is.null(x = ncol)) {
              ncol <- 2
              if (length(x = features) == 1) {
                  ncol <- 1
              }
              if (length(x = features) > 6) {
                  ncol <- 3
              }
              if (length(x = features) > 9) {
                  ncol <- 4
              }
          }
          ncol <- ifelse(test = is.null(x = split.by) || blend,
              yes = ncol, no = length(x = features))
          legend <- if (blend) {
              "none"
          }
          else {
              split.by %iff% "none"
          }
          plots <- CombinePlots(plots = plots, ncol = ncol, legend = legend,
              nrow = split.by %iff% length(x = levels(x = data$split)))
      }
      return(plots)
}
environment(XY_FeaturePlot_no_order) <- asNamespace('Seurat')

runFDG = function(pca.df, snn, iterations = 600, tool_addr, python.addr){
  current.wd = getwd()
  # generate unique name for pca data file
  pca.data.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  pca.data.fname = paste(pca.data.fname, ".csv", sep = "")
  # generate unique name for snn file
  snn.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  snn.fname = paste(snn.fname, ".smm", sep = "")
  # generate unique name for fdg coordinates
  fdg.coordinates.fname = paste(sample(LETTERS, 20, TRUE), collapse = "")
  fdg.coordinates.fname = paste(fdg.coordinates.fname, ".csv", sep = "")
  write.csv(pca.df, pca.data.fname)
  writeMM(obj=snn, file=snn.fname)
  command = gsub(pattern="ITER", replacement=as.character(iterations), paste(tool_addr, "ITER", sep = " "))
  command = paste(command, paste(c(pca.data.fname, snn.fname, fdg.coordinates.fname), collapse = " "), sep = " ")
  system(command, wait = T)
  fdg_coordinates = read.csv(fdg.coordinates.fname, header = FALSE)
  colnames(fdg_coordinates) = c("X", "Y")
  rownames(fdg_coordinates) = rownames(pca.df)
  file.remove(c(pca.data.fname, snn.fname, fdg.coordinates.fname))
  setwd(current.wd)
  return(fdg_coordinates)
}

seuratToURD2 <- function(seurat.object) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    
    # Copy over data
    ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
    if(!any(dim(seurat.object@assays$RNA@counts) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$RNA@counts[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
    # Copy over metadata
    ## TO DO - grab kmeans clustering info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) { 
      get.data <- as.data.frame(seurat.object@assays$RNA@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) { 
      get.data <- as.data.frame(seurat.object@meta.data) 
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    # Copy over var.genes
    if(length(seurat.object@assays$RNA@var.features > 0)) ds@var.genes <- seurat.object@assays$RNA@var.features
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    }
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("pca" %in% names(seurat.object@reductions)) && !any(dim(Loadings(seurat.object, reduction = "pca")) == 0)) {
        ds@pca.load <- as.data.frame(Loadings(seurat.object, reduction = "pca"))
        ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
        ds@pca.sdev <- seurat.object@reductions$pca@stdev
        ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
}
environment(seuratToURD2) <- asNamespace('URD')

XY_PlotUTRLengthShift <- function (results.table, plot.title = "Global shift in 3'UTR length",
    do.ranksum.test = TRUE, return.plot = TRUE, do.plot = FALSE,binwidth,base_size)
{
    locations.results.table.up <- subset(results.table, FC_direction == "Up")
    pos.upreg <- apply(as.matrix(locations.results.table.up[, c("SiteLocation",
        "NumSites")]), 1, function(x) {
        relative_location(x[1], x[2])
    })
    locations.results.table.down <- subset(results.table, FC_direction ==
        "Down")
    pos.downreg <- apply(as.matrix(locations.results.table.down[,
        c("SiteLocation", "NumSites")]), 1, function(x) {
        relative_location(x[1], x[2])
    })
    if (do.ranksum.test) {
        this.test <- wilcox.test(pos.upreg, pos.downreg)
        print("Wilcoxon Rank-sum test comparing relative peak locations for up- vs down-regulated peaks:")
        print(paste0("P-value = ", this.test$p.value))
    }
    ggData <- data.frame(Peak_location = c(pos.upreg, pos.downreg),
        FC_direction = c(rep("Up", length(pos.upreg)), rep("Down",
            length(pos.downreg))))
    ggData$FC_direction <- factor(ggData$FC_direction, levels = c("Up",
        "Down"))
    pl.density <- ggplot(ggData, aes(Peak_location, stat(count),
        fill = FC_direction)) + geom_density(alpha = 0.8) + ylab("") +
        theme_void(base_size = base_size) + theme(axis.text = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank()) +
        ggtitle(plot.title) + theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette = "Set1")
    pl.histogram <- ggplot(ggData, aes(Peak_location, fill = FC_direction)) +
        geom_histogram(position = position_dodge(), colour = "black",
            binwidth = binwidth, alpha = 0.8) + theme_classic(base_size = base_size) +
        xlab("Relative peak location") + ylab("Peak count") +
        scale_fill_brewer(palette = "Set1") + theme(legend.position = "right") +
        guides(fill = guide_legend(title = "Fold-change\ndirection",
            title.position = "top"))
    pl.combined <- cowplot::plot_grid(pl.density, pl.histogram,
        ncol = 1, rel_heights = c(0.3, 0.8), axis = "lr", align = "v")
    if (do.plot) {
        plot(pl.combined)
    }
    if (return.plot) {
        return(pl.combined)
    }
}
environment(XY_PlotUTRLengthShift) <- asNamespace("Sierra")

XY_heatmap <- function (seurat_obj=seurat_obj, group=group,genes=genes,all_num=all_num,assay_sel=assay_sel,labels_rot=labels_rot,
  color=color,min_and_max_cut=num_cut,new_names=new_names,show_row_names=show_row_names,mark_gene=mark_gene,label_size=label_size,scale=scale){
  message("Processed data begain")
  ATAC <- GetAssayData(seurat_obj,slot="data",assay=assay_sel)
  ATAC_sel <- ATAC[genes,]
  ATAC_sel <- as.matrix(ATAC_sel)
  if (scale==TRUE) {
    ATAC_sel_zscore <- t(apply(ATAC_sel, 1, function(x) (x-mean(x))/sd(x)))
    } else {
      ATAC_sel_zscore <- ATAC_sel
    }
  sel_cutoff <- min(abs(range(ATAC_sel_zscore)))
  if (is.null(min_and_max_cut)){
    ATAC_sel_zscore[ATAC_sel_zscore > sel_cutoff] <- sel_cutoff
    ATAC_sel_zscore[ATAC_sel_zscore < -sel_cutoff] <- -sel_cutoff
    } else {
      ATAC_sel_zscore[ATAC_sel_zscore > min_and_max_cut] <- min_and_max_cut
      ATAC_sel_zscore[ATAC_sel_zscore < -min_and_max_cut] <- -min_and_max_cut
    }
  meta_info <- seurat_obj@meta.data
  meta_info <- meta_info[order(meta_info[,group],decreasing=F),]
  annotation = data.frame(new_anno=meta_info[group],cell_id=rownames(meta_info))
  colnames(annotation) <- c("new_anno","cell_id")
  if (all_num == TRUE & !is.null(new_names)){
    annotation$new_anno <- paste(new_names,annotation$new_anno,sep="")
    aa <- as.data.frame(table(annotation$new_anno))
    aa <- aa[order(aa$Freq,decreasing=T),]
    aa$Var1 <- as.character(aa$Var1)
    annotation$new_anno <- factor(annotation$new_anno,levels=aa$Var1)
  }
  annotation <- annotation[order(annotation$new_anno),]
  annotation = data.frame(new_anno=annotation$new_anno,row.names=rownames(annotation))
  require(pheatmap)
  message("pheatmap printing start")
  require(ComplexHeatmap)
  require(BuenColors)
  require(scales) 
  col1 <- jdb_palette("Darjeeling2")
  col2 <- jdb_palette("Darjeeling")
  col3 <- jdb_palette("Moonrise3")
  col_sel <- c(col1,col2,col3)
  col_sel <- hue_pal()(length(as.character(unique(annotation$new_anno))))
  col <- col_sel[1:length(as.character(unique(annotation$new_anno)))]
  names(col) <- as.character(unique(annotation$new_anno))
  top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # 设置填充色
  labels = as.character(unique(annotation$new_anno)), 
  labels_gp = gpar(cex = label_size , col = "black"),labels_rot=labels_rot))
  if (is.null(mark_gene)){
  ph <- Heatmap(ATAC_sel_zscore[,rownames(annotation)],
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = show_row_names,
    top_annotation = top_anno,
    col = rev(color),
    column_split = annotation$new_anno,
    column_title_rot = 90)
    } else {
      both_gene <- intersect(rownames(ATAC_sel_zscore[,rownames(annotation)]),mark_gene)
      gene_pos <- which(rownames(ATAC_sel_zscore[,rownames(annotation)]) %in% both_gene)
      selected_gene <- rownames(ATAC_sel_zscore[,rownames(annotation)])[gene_pos]
      row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = selected_gene))
      ph <- Heatmap(ATAC_sel_zscore[,rownames(annotation)],
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = show_row_names,
        top_annotation = top_anno,
        right_annotation = row_anno,
        col = rev(color),
        column_split = annotation$new_anno,
        column_title_rot = 90)
    }
    return(ph)
    print(ph)
}


anno_scRNA <- function(refer_pos,number.group,outpath,outnames){
  reference_tmp <- list.files(refer_pos)
  reference_files <- paste0(refer_pos,reference_tmp)
  names <- gsub(".csv","",reference_tmp)
  require(dplyr)
  score_predict_all_tmp <- lapply(1:length(reference_files),function(which_file){
    reference_files_selec <- as.character(reference_files[which_file])
    reference_files_selec <- read.csv(file=reference_files_selec, header=T)
    reference_files_selec <- arrange(reference_files_selec, cluster, desc(p_val))
    cell_types <- unique(as.character(reference_files_selec$Annotation))
    all_score_tmp <- lapply(cell_types,function(i){
      y <- reference_files_selec$Annotation
      specific_symbols <- subset(reference_files_selec,Annotation==i)
      specific_symbols$order <- c(1:nrow(specific_symbols))
      specific_score <- lapply(as.numeric(unique(as.character(number.group))),function(query_group){
        cluster <- subset(top200,cluster==query_group)
        cluster_specific_symbols <- merge(cluster,specific_symbols, by = "gene")
        omit_cluster_specific <- na.omit(cluster_specific_symbols)
        cluster_specific_score  <- (200*sum(omit_cluster_specific$order))/(nrow(specific_symbols)*(nrow(specific_symbols)+1))
        cluster_specific_score <- as.data.frame(cluster_specific_score)
        rownames(cluster_specific_score) <- paste("cluster",query_group,sep="_")
        return(cluster_specific_score)
        })
      specific_score_all_cluster <- do.call(rbind,specific_score)
      colnames(specific_score_all_cluster) <- paste0(i,"score")
      specific_score_all_cluster <- t(specific_score_all_cluster)
      return(specific_score_all_cluster)
      })
    all_score <- do.call(rbind,all_score_tmp)
    all_score <- round(all_score,3)
    all_score <- as.data.frame(t(all_score))
    all_score_predi <- lapply(1:nrow(all_score),function(rows){
      sel_rows <- all_score[rows,]
      max_names <- colnames(sel_rows)[which(sel_rows==max(sel_rows))]
      sel_rows$MAX_Predic <- max_names
      sel_rows$MAX_Score <- sel_rows[,max_names]
      return(sel_rows)
      })
    all_score_predi <- do.call(rbind,all_score_predi)
    all_score_predi$referen_orig <- names[which_file]
    hh <- paste("calculated_Score_from",which_file,".csv",sep="_")
    if (length(outnames)==0){
      outpath <- outpath
      cc <- paste(outpath,hh,sep="/")
    } else {
      outpath <- outpath
      cc <- paste(outpath,outnames,hh,sep="/")
    }
    write.csv(all_score_predi,file=cc)
    all_score_ext <- data.frame(Query_Clu=rownames(all_score_predi),all_score_predi$MAX_Predic,all_score_predi$MAX_Score,all_score_predi$referen_orig)
    message(names[which_file]," is done")
    return(all_score_ext)
    })
  score_predict_all <- do.call(rbind,score_predict_all_tmp)
  if (length(outnames)==0){
    outpath <- outpath
    cc <- paste(outpath,"_Summery.csv",sep="/")
  } else {
    outpath <- outpath
    cc <- paste(outpath,outnames,"_Summery.csv",sep="/")
  }
  write.csv(score_predict_all,file=cc)
  message("Summery is done")
}



XY_simulateRandomWalksFromTips <- function (object, tip.group.id, root.cells, transition.matrix,
    n.per.tip = 10000, root.visits = 1, max.steps = ncol(object@logupx.data),
    verbose = T)
{
  require(future)
  require(future.apply)
  options(future.globals.maxSize = 3000 * 1024^2)
  plan("multiprocess", workers = 8)
  plan()
    all.tips <- sort(setdiff(as.character(unique(object@group.ids[,
        tip.group.id])), NA))
    n.tips <- length(all.tips)
    tip.walks <- future_lapply(all.tips, function(tip) {
        if (verbose)
            print(paste0(Sys.time(), " - Starting random walks from tip ",
                tip))
        tip.cells <- rownames(object@group.ids)[which(as.character(object@group.ids[,
            tip.group.id]) == tip)]
        if (verbose)
            verbose.freq = round(n.per.tip/10)
        else verbose.freq = 0
        walks <- simulateRandomWalk(start.cells = tip.cells,
            transition.matrix = transition.matrix, end.cells = root.cells,
            n = n.per.tip, end.visits = root.visits, verbose.freq = verbose.freq,
            max.steps = max.steps)
        return(walks)
    })
    names(tip.walks) <- all.tips
    return(tip.walks)
}
environment(XY_simulateRandomWalksFromTips) <- asNamespace('URD')



XY_floodPseudotimeProcess <- function (object, floods, floods.name = "pseudotime", max.frac.NA = 0.4,
    pseudotime.fun = mean, stability.div = 10)
{
    if (nrow(object@diff.data) == 0) {
        warning("Initializing @diff.data, though this should have been previously initialized by creation or importation of diffusion map.")
        if (nrow(object@dm@eigenvectors) == 0) {
            stop("Make sure your diffusion map is loaded into the object (calcDM or importDM)")
        }
        object@diff.data <- data.frame(row.names = rownames(dm@eigenvectors))
    }
    if (nrow(object@pseudotime) == 0) {
        warning("Initializing @pseudotime, though this should have been previously initialized by creation or importation of diffusion map.")
        if (nrow(object@dm@eigenvectors) == 0) {
            stop("Make sure your diffusion map is loaded into the object (calcDM or importDM)")
        }
        object@pseudotime <- data.frame(row.names = rownames(dm@eigenvectors))
    }
    floods.name <- as.character(floods.name)
    if (class(floods) == "list") {
        floods <- do.call("cbind", floods)
    }
    else if (class(floods) != "data.frame") {
        stop("floods should be a list or data.frame.")
    }
    require(future)
    require(future.apply)
    options(future.globals.maxSize = 3000 * 1024^2)
    plan("multiprocess", workers = 8)
    plan()
    frac.na <- future_apply(floods, 1, function(x) sum(is.na(x)))/dim(floods)[2]
    cells.toss <- names(which(frac.na > max.frac.NA))
    floods <- floods[setdiff(rownames(floods), cells.toss), ]
    floods <- sweep(floods, 2, apply(floods, 2, max, na.rm = T),
        "/")
    floods.in.division <- ceiling(1:stability.div/stability.div *
        dim(floods)[2])
    walks.per.cell <- as.data.frame(lapply(floods.in.division,
        function(n.floods) {
            if (n.floods == 1)
                return(as.numeric(!is.na(floods[, 1])))
            else return(apply(floods[, 1:n.floods], 1, function(x) sum(as.numeric(!is.na(x)))))
        }))
    pseudotime.stability <- as.data.frame(lapply(floods.in.division,
        function(n.floods) {
            if (n.floods == 1)
                return(floods[, 1])
            else return(apply(floods[, 1:n.floods], 1, mean,
                na.rm = T))
        }))
    colnames(walks.per.cell) <- seq(1, ncol(floods), length.out = stability.div)
    colnames(pseudotime.stability) <- seq(1, ncol(floods), length.out = stability.div)
    object@pseudotime.stability$pseudotime <- pseudotime.stability
    object@pseudotime.stability$walks.per.cell <- walks.per.cell
    final.visit.freq <- walks.per.cell[rownames(object@diff.data),
        dim(walks.per.cell)[2]]
    object@diff.data[, paste0("visitfreq.raw.", floods.name)] <- final.visit.freq
    object@diff.data[, paste0("visitfreq.log.", floods.name)] <- log10(final.visit.freq +
        1)
    object@pseudotime[, floods.name] <- NA
    object@pseudotime[rownames(pseudotime.stability), floods.name] <- pseudotime.stability[,
        stability.div]
    return(object)
}
environment(XY_simulateRandomWalksFromTips) <- asNamespace('URD')



XY_treeForceDirectedLayout <- function (object, num.nn = NULL, method = c("fr", "drl", "kk"),
    cells.to.do = NULL, cell.minimum.walks = 1, cut.outlier.cells = NULL,
    cut.outlier.edges = NULL, max.pseudotime.diff = NULL, cut.unconnected.segments = 2,
    min.final.neighbors = 2, remove.duplicate.cells = T, tips = object@tree$tips,
    coords = "auto", start.temp = NULL, n.iter = NULL, density.neighbors = 10,
    plot.outlier.cuts = F, verbose = F)
{
    require(future)
    require(future.apply)
    options(future.globals.maxSize = 3000 * 1024^2)
    plan("multiprocess", workers = 8)
    plan()
    if (length(method) > 1)
        method <- method[1]
    if (is.null(cells.to.do))
        cells.to.do <- rownames(object@diff.data)
    starting.cells <- length(cells.to.do)
    if (is.null(num.nn))
        num.nn <- ceiling(sqrt(length(cells.to.do)))
    if (is.null(tips))
        tips <- object@tree$tips
    if (class(coords) == "character" && coords == "auto") {
        coords <- as.matrix(object@tree$cell.layout[cells.to.do,
            c("x", "y")])
    }
    else if (!is.null(coords) && class(coords) != "matrix") {
        stop("coords must either be 'auto', NULL, or a matrix.")
    }
    dim = 2
    if (verbose)
        print(paste(Sys.time(), ": Starting with parameters",
            method, num.nn, "NN", dim, "D", length(cells.to.do),
            "cells"))
    cells.no.pseudotime <- cells.to.do[which(is.na(object@tree$pseudotime[cells.to.do]))]
    if (!is.na(cut.unconnected.segments) & !is.null(cut.unconnected.segments) &
        cut.unconnected.segments > 0) {
        cells.not.in.tree <- setdiff(cells.to.do, unlist(object@tree$cells.in.segment))
        if (verbose)
            print(paste0("Removing ", length(unique(c(cells.no.pseudotime,
                cells.not.in.tree))), " cells that are not assigned a pseudotime or a segment in the tree."))
        cells.to.do <- setdiff(cells.to.do, c(cells.no.pseudotime,
            cells.not.in.tree))
    }
    else {
        if (verbose)
            print(paste0("Removing ", length(cells.no.pseudotime),
                " cells that are not assigned a pseudotime."))
        cells.to.do <- setdiff(cells.to.do, cells.no.pseudotime)
    }
    if (verbose)
        print(paste0(Sys.time(), ": Preparing walk data."))
    walk.data <- object@diff.data[cells.to.do, paste0("visitfreq.raw.",
        object@tree$tips)]
    walk.data$pseudotime <- object@tree$pseudotime[cells.to.do]
    walk.total <- future_apply(walk.data, 1, sum)
    if (any(walk.total < cell.minimum.walks)) {
        if (verbose)
            print(paste0("Removing ", length(which(walk.total <
                cell.minimum.walks)), " cells that were visited fewer than ",
                cell.minimum.walks, " times by random walks."))
        walk.data <- walk.data[which(walk.total >= cell.minimum.walks),
            ]
        walk.total <- future_apply(walk.data, 1, sum)
    }
    walk.data <- sweep(walk.data, 1, walk.total, "/")
    walk.data <- as.matrix(walk.data)
    duped <- which(duplicated(walk.data))
    if (length(duped) > 0 && remove.duplicate.cells) {
        warning(length(duped), " cells have duplicate random walk coordinates and are being removed from the layout.")
        walk.data <- unique(walk.data)
    }
    else if (length(duped) > 0) {
        warning(length(duped), " cells have duplicate random walk coordinates. This may cause a problem in the layout. If so, set remove.duplicate.cells=T")
    }
    if (verbose)
        print(paste0(Sys.time(), ": Calculating nearest neighbor graph."))
    walk.nn <- RANN::nn2(data = walk.data, query = walk.data,
        k = max(num.nn) + 1, treetype = "kd", searchtype = "priority")
    rownames(walk.nn$nn.idx) <- rownames(walk.data)
    walk.nn$nn.idx <- walk.nn$nn.idx[, 2:(num.nn + 1)]
    rownames(walk.nn$nn.dists) <- rownames(walk.data)
    walk.nn$nn.dists <- walk.nn$nn.dists[, 2:(num.nn + 1)]
    walk.nn$nn.label <- future_apply(walk.nn$nn.idx, 2, function(y) rownames(walk.data)[y])
    rownames(walk.nn$nn.label) <- rownames(walk.data)
    if (!is.null(cut.outlier.cells)) {
        q <- quantile(walk.nn$nn.dists[, 2])
        outer.fence <- q[4] + cut.outlier.cells * (q[4] - q[1])
        if (plot.outlier.cuts) {
            hist(walk.nn$nn.dists[, 2], breaks = 100, main = "Outlier Removal: Second Nearest Neighbor",
                xlab = "Distance to 2nd Nearest Neighbor")
            abline(v = outer.fence, col = "red")
        }
        if (verbose)
            print(paste0(Sys.time(), ": Removing ", round(length(which(walk.nn$nn.dists[,
                2] > outer.fence))/starting.cells, digits = 2),
                "% of cells as outliers."))
        cells.keep <- which(walk.nn$nn.dists[, 2] <= outer.fence)
    }
    else {
        cells.keep <- 1:dim(walk.data)[1]
    }
    if (verbose)
        print(paste0(Sys.time(), ": Preparing edge list."))
    edges <- reshape2::melt(walk.nn$nn.label[cells.keep, ], stringsAsFactors = F)
    dists <- reshape2::melt(walk.nn$nn.dists[cells.keep, ])
    edges <- edges[, c(1, 3)]
    names(edges) <- c("V1", "V2")
    edges$dists <- dists$value
    edges$V1 <- as.character(edges$V1)
    edges$V2 <- as.character(edges$V2)
    starting.edges <- dim(edges)[1]
    if (!is.null(cut.outlier.edges)) {
        q <- quantile(edges$dists)
        outer.fence <- q[4] + cut.outlier.edges * (q[4] - q[1])
        if (verbose)
            print(paste0(Sys.time(), ": Removing ", round(length(which(edges$dists >
                outer.fence))/starting.edges * 100, digits = 3),
                "% of edges as outliers."))
        if (plot.outlier.cuts) {
            hist(edges$dists, breaks = 100, main = "Outlier Removal: Edge distances",
                xlab = "Edge distance")
            abline(v = outer.fence, col = "red")
        }
        edges <- edges[edges$dists <= outer.fence, ]
    }
    if (!is.null(max.pseudotime.diff)) {
        edges$pt1 <- object@tree$pseudotime[edges$V1]
        edges$pt2 <- object@tree$pseudotime[edges$V2]
        edges$dpt <- abs(edges$pt1 - edges$pt2)
        if (verbose)
            print(paste0(Sys.time(), ": Removing ", round(length(which(edges$dpt >
                max.pseudotime.diff))/starting.edges * 100, digits = 3),
                "% of edges with too large pseudotime difference."))
        edges <- edges[which(edges$dpt <= max.pseudotime.diff),
            ]
    }
    if (!is.na(cut.unconnected.segments) & !is.null(cut.unconnected.segments) &
        cut.unconnected.segments > 0) {
        edges$seg1 <- object@diff.data[edges$V1, "segment"]
        edges$seg2 <- object@diff.data[edges$V2, "segment"]
        seg.dist <- igraph::distances(object@tree$tree.igraph,
            mode = "all")
        edges$seg.dist <- future_apply(edges[, c("seg1", "seg2")], 1,
            function(segs) seg.dist[segs[1], segs[2]])
        if (verbose)
            print(paste0(Sys.time(), ": Removing ", round(length(which(edges$seg.dist >
                cut.unconnected.segments))/starting.edges * 100,
                digits = 2), "% of edges that are between segments with distance > ",
                cut.unconnected.segments))
        edges <- edges[edges$seg.dist <= cut.unconnected.segments,
            ]
    }
    if (verbose)
        print(paste0(Sys.time(), ": Trimming cells that are no longer well connected."))
    connections.remaining <- table(edges[edges$dists > 0, "V1"])
    cells.with.enough.connections <- names(which(connections.remaining >=
        min.final.neighbors))
    cells.without.enough.connections <- length(unique(edges$V1)) -
        length(cells.with.enough.connections)
    while (cells.without.enough.connections > 0) {
        edges <- edges[which(edges$V1 %in% cells.with.enough.connections &
            edges$V2 %in% cells.with.enough.connections), ]
        connections.remaining <- table(edges[edges$dists > 0,
            "V1"])
        cells.with.enough.connections <- names(which(connections.remaining >=
            min.final.neighbors))
        cells.without.enough.connections <- length(unique(edges$V1)) -
            length(cells.with.enough.connections)
    }
    if (verbose)
        print(paste0(Sys.time(), ": ", round(length(cells.with.enough.connections)/starting.cells *
            100, digits = 2), "% of starting cells preserved."))
    if (verbose)
        print(paste0(Sys.time(), ": Preparing igraph object."))
    if (method == "kk") {
        edges$weight <- edges$dists
    }
    else {
        edges$weight <- 1/edges$dists
    }
    object@tree$walks.force.edges <- edges
    igraph.walk.weights <- igraph::graph_from_data_frame(edges,
        directed = F)
    if (!is.null(coords)) {
        cells.remain <- unique(c(edges$V1, edges$V2))
        coords.remain <- rownames(coords) %in% cells.remain
        coords <- coords[coords.remain, ]
    }
    if (method == "drl" && !igraph::is_connected(igraph.walk.weights)) {
        igraph.walk.weights.decomposed <- decompose(igraph.walk.weights)
        graph.keep <- which.max(unlist(future_lapply(igraph.walk.weights.decomposed,
            function(x) length(V(x)))))
        warning("DrL method can only operate on fully connected graphs. Using more nearest neighbors will increase graph connectivity. For now, using largest subgraph: ",
            length(V(igraph.walk.weights.decomposed[[graph.keep]])),
            " of ", length(V(igraph.walk.weights)), " cells.")
        igraph.walk.weights <- igraph.walk.weights.decomposed[[graph.keep]]
        coords <- coords[rownames(coords) %in% names(V(igraph.walk.weights)),
            ]
        cells.with.enough.connections <- intersect(cells.with.enough.connections,
            names(V(igraph.walk.weights)))
    }
    if (is.null(start.temp))
        start.temp <- sqrt(igraph::vcount(igraph.walk.weights))
    if (is.null(n.iter)) {
        if (method == "fr")
            n.iter = 500
        if (method == "kk")
            n.iter = 50 * igraph::vcount(igraph.walk.weights)
    }
    if (verbose)
        print(paste0(Sys.time(), ": Doing force-directed layout."))
    if (method == "fr")
        igraph.walk.layout <- igraph::layout_with_fr(igraph.walk.weights,
            dim = dim, coords = coords, start.temp = start.temp,
            niter = n.iter)
    if (method == "drl")
        igraph.walk.layout <- igraph::layout_with_drl(igraph.walk.weights,
            dim = dim, options = list(edge.cut = 0))
    if (method == "kk")
        igraph.walk.layout <- igraph::layout_with_kk(igraph.walk.weights,
            dim = dim, coords = coords, maxiter = n.iter)
    if (dim == 2) {
        if (verbose)
            print(paste0(Sys.time(), ": Calculating Z."))
        neighbor.pt <- unlist(future_lapply(cells.with.enough.connections,
            function(i) {
                weights <- 1/walk.nn$nn.dists[i, ]
                weights[!is.finite(weights)] <- sum(weights[is.finite(weights)])
                weighted.mean(x = object@tree$pseudotime[walk.nn$nn.label[i,
                  ]], w = weights)
            }))
        names(neighbor.pt) <- cells.with.enough.connections
        object@tree$walks.force.layout <- as.data.frame(igraph.walk.layout,
            stringsAsFactors = F)
        names(object@tree$walks.force.layout) <- c("x", "y")
        rownames(object@tree$walks.force.layout) <- igraph::V(igraph.walk.weights)$name
        fdl.failed <- length(which(!(complete.cases(object@tree$walks.force.layout))))
        if (fdl.failed > 0) {
            warning(fdl.failed, " cells were not assigned coordinates during the force-directed layout.")
            object@tree$walks.force.layout <- object@tree$walks.force.layout[complete.cases(object@tree$walks.force.layout),
                ]
        }
        neighbor.pt.factor <- mean(future_apply(object@tree$walks.force.layout[,
            c("x", "y")], 2, max))/max(neighbor.pt, na.rm = T)
        object@tree$walks.force.layout$telescope.pt <- neighbor.pt[rownames(object@tree$walks.force.layout)] *
            neighbor.pt.factor
        if (verbose)
            print(paste0(Sys.time(), ": Calculating local density."))
        object <- fdlDensity(object, neighbor = density.neighbors)
    }
    else if (dim == 3) {
        object@tree$walks.force.layout.3d <- as.data.frame(igraph.walk.layout,
            stringsAsFactors = F)
        names(object@tree$walks.force.layout.3d) <- c("x", "y",
            "z")
        rownames(object@tree$walks.force.layout.3d) <- igraph::V(igraph.walk.weights)$name
    }
    if (!is.null(object@tree$segment.names)) {
        object@tree$walks.force.labels <- treeForcePositionLabels(object)
    }
    if (verbose)
        print(paste0(Sys.time(), ": Finished."))
    return(object)
}
environment(XY_treeForceDirectedLayout) <- asNamespace('URD')



XY_graphClustering <- function (object, dim.use = c("pca", "dm"), cells.use = NULL,
    which.dims = which(object@pca.sig), num.nn = 30, do.jaccard = TRUE,
    method = c("Louvain", "Infomap"), group.id = method)
{
    require(future)
    require(future.apply)
    options(future.globals.maxSize = 3000 * 1024^2)
    plan("multiprocess", workers = 8)
    plan()
    if (length(dim.use) > 1)
        dim.use <- dim.use[1]
    if (length(method) > 1)
        method <- method[1]
    if (dim.use == "pca") {
        if (is.null(cells.use)) {
            data.use = object@pca.scores[, which.dims]
        }
        else {
            data.use = object@pca.scores[cells.use, which.dims]
        }
    }
    else if (dim.use == "dm") {
        if (is.null(cells.use)) {
            data.use = object@dm@eigenvectors[, which.dims]
        }
        else {
            data.use = object@dm@eigenvectors[cells.use, which.dims]
        }
    }
    else {
        stop("dim.use must be either pca or dm")
    }
    nearest <- nn2(data.use, data.use, k = max(num.nn) + 1, treetype = "bd",
        searchtype = "priority")
    nearest$nn.idx <- nearest$nn.idx[, -1]
    nearest$nn.dists <- nearest$nn.dists[, -1]
    for (this.nn in num.nn) {
        edges = melt(t(nearest$nn.idx[, 1:this.nn]))
        colnames(edges) = c("B", "A", "C")
        edges = edges[, c("A", "B", "C")]
        edges$B = edges$C
        edges$C = 1
        edges = unique(transform(edges, A = pmin(A, B), B = pmax(A,
            B)))
        if (do.jaccard) {
            NN = nearest$nn.idx[, 1:this.nn]
            jaccard_dist = future_apply(edges, 1, function(x) length(intersect(NN[x[1],
                ], NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2],
                ])))
            edges$C = jaccard_dist
            edges = subset(edges, C != 0)
            edges$C = edges$C/max(edges$C)
        }
        names(edges) <- c("V1", "V2", "weight")
        edges$V1 <- rownames(data.use)[edges$V1]
        edges$V2 <- rownames(data.use)[edges$V2]
        g <- graph.data.frame(edges, directed = F)
        if (tolower(method) == "louvain") {
            graph.out = cluster_louvain(g)
        }
        else if (tolower(method) == "infomap") {
            graph.out = cluster_infomap(g)
        }
        else {
            stop("Method should be either Louvain or Infomap.")
        }
        clust.assign = factor(graph.out$membership, levels = sort(unique(graph.out$membership)))
        names(clust.assign) = graph.out$names
        k = order(table(clust.assign), decreasing = TRUE)
        new.levels = rep(1, length(unique(graph.out$membership)))
        new.levels[k] = 1:length(unique(graph.out$membership))
        levels(clust.assign) = new.levels
        clust.assign = factor(clust.assign, levels = 1:length(unique(graph.out$membership)))
        object@meta$clust = NULL
        object@meta[names(clust.assign), "clust"] = clust.assign
        this.group.id <- paste0(group.id, "-", this.nn)
        object@group.ids[names(clust.assign), this.group.id] <- as.character(clust.assign)
    }
    return(object)
}
environment(XY_graphClustering) <- asNamespace('URD')


XY_processRandomWalks <- function (object, walks, walks.name, aggregate.fun = mean, n.subsample = 10,verbose = T)
{
    require(future)
    require(future.apply)
    options(future.globals.maxSize = 3000 * 1024^2)
    plan("multiprocess", workers = 8)
    plan()
    walks.name <- as.character(walks.name)
    walks <- walks[!unlist(future_lapply(walks, is.null))]
    hops.melt <- as.data.frame(data.table::rbindlist(future_lapply(walks,
        function(i) {
            data.frame(hops = (1:length(i))/length(i), cell = i,
                stringsAsFactors = F)
        })), stringsAsFactors = F)
    walks.in.division <- ceiling(1:n.subsample/n.subsample *
        length(walks))
    walk.lengths <- unlist(future_lapply(walks, length))
    pseudotime.stability <- matrix(rep(NA, dim(object@diff.data)[1] *
        n.subsample), nrow = dim(object@diff.data)[1], ncol = n.subsample,
        dimnames = list(rownames(object@diff.data), walks.in.division))
    walks.per.cell <- matrix(rep(0, dim(object@diff.data)[1] *
        n.subsample), nrow = dim(object@diff.data)[1], ncol = n.subsample,
        dimnames = list(rownames(object@diff.data), walks.in.division))
    for (section in 1:n.subsample) {
        if (verbose)
            print(paste("Calculating pseudotime with", walks.in.division[section],
                "walks."))
        visit.freq <- table(unlist(walks[1:walks.in.division[section]]))
        walks.per.cell[, section] <- visit.freq[rownames(walks.per.cell)]
        hops.melt.rows <- sum(walk.lengths[1:walks.in.division[section]])
        hops.relative <- reshape2::dcast(data = hops.melt[1:hops.melt.rows,
            ], formula = cell ~ ., fun.aggregate = aggregate.fun,
            value.var = "hops")
        rownames(hops.relative) <- hops.relative$cell
        pseudotime.stability[, section] <- hops.relative[rownames(pseudotime.stability),
            "."]
    }
    object@pseudotime.stability$pseudotime <- pseudotime.stability
    object@pseudotime.stability$walks.per.cell <- walks.per.cell
    final.visit.freq <- walks.per.cell[rownames(object@diff.data),
        dim(walks.per.cell)[2]]
    final.visit.freq[is.na(final.visit.freq)] <- 0
    object@diff.data[, paste0("visitfreq.raw.", walks.name)] <- final.visit.freq
    object@diff.data[, paste0("visitfreq.log.", walks.name)] <- log10(final.visit.freq +
        1)
    object@pseudotime[, walks.name] <- NA
    object@pseudotime[hops.relative[, "cell"], walks.name] <- hops.relative[,
        "."]
    return(object)
}
environment(XY_processRandomWalks) <- asNamespace('URD')

XY_processRandomWalksFromTips <- function (object, walks.list, aggregate.fun = mean, n.subsample = 10,verbose = T)
{
    if (class(walks.list) != "list" || class(walks.list[[1]]) !=
        "list")
        stop("walks.list must be a list (tips) of lists (walks).")
    if (is.null(names(walks.list)))
        stop("walks.list must be named according to tips.")
    tip.names <- names(walks.list)
    for (tip in tip.names) {
        if (verbose)
            print(paste0(Sys.time(), " - Processing walks from tip ",
                tip))
        object <- XY_processRandomWalks(object, walks = walks.list[[tip]],
            walks.name = tip, aggregate.fun = aggregate.fun,
            n.subsample = n.subsample, verbose = verbose)
    }
    return(object)
}
environment(XY_processRandomWalksFromTips) <- asNamespace('URD')




XY_buildTree <- function (object, pseudotime, tips.use = NULL, divergence.method = c("ks",
    "preference"), weighted.fusion = T, use.only.original.tips = T,
    cells.per.pseudotime.bin = 80, bins.per.pseudotime.window = 5,
    minimum.visits = 10, visit.threshold = 0.7, save.breakpoint.plots = NULL,
    save.all.breakpoint.info = F, p.thresh = 0.01, min.cells.per.segment = 1,
    min.pseudotime.per.segment = 0.01, dendro.node.size = 100,
    dendro.cell.jitter = 0.15, dendro.cell.dist.to.tree = 0.05,
    verbose = T)
{
    require(future)
    require(future.apply)
    options(future.globals.maxSize = 3000 * 1024^2)
    plan("multiprocess", workers = 8)
    plan()
    if (length(divergence.method) > 1)
        divergence.method <- divergence.method[1]
    if (!(divergence.method %in% c("ks", "preference")))
        stop("Divergence method must be 'ks' or 'preference'.")
    tips <- as.character(tips.use)
    object@tree$tips <- tips
    if (weighted.fusion) {
        tip.size <- unlist(future_lapply(object@tree$cells.in.tip, length))
        if (length(tip.size) == 0)
            stop("Either weighted.fusion must be false, or loadTipCells must be run prior.")
    }
    object@tree$pseudotime <- object@pseudotime[, pseudotime]
    names(object@tree$pseudotime) <- rownames(object@pseudotime)
    seg.add <- as.character(max(suppressWarnings(as.numeric(tips))) +
        1)
    if (is.na(seg.add))
        seg.add <- "1"
    cells.in.segments <- putativeCellsInSegment(object, tips,
        minimum.visits = minimum.visits, visit.threshold = visit.threshold)
    object@tree$cells.in.segments <- cells.in.segments
    object@tree$segment.pseudotime.limits <- data.frame(start = future_sapply(tips,
        function(segment) min(object@pseudotime[rownames(cells.in.segments)[which(cells.in.segments[,
            segment])], pseudotime])), end = future_sapply(tips, function(segment) max(object@pseudotime[rownames(cells.in.segments)[which(cells.in.segments[,
        segment])], pseudotime])), stringsAsFactors = F, row.names = tips)
    object <- allSegmentDivergenceByPseudotime(object, pseudotime = pseudotime,
        divergence.method = divergence.method, segments = tips,
        pseudotime.cuts = cells.per.pseudotime.bin, window.size = bins.per.pseudotime.window,
        minimum.visits = minimum.visits, visit.threshold = visit.threshold,
        p.thresh = p.thresh, breakpoint.decision.plots = save.breakpoint.plots,
        cache = F, verbose = verbose)
    if (save.all.breakpoint.info) {
        all.pseudotime.breakpoint.details <- object@tree$pseudotime.breakpoint.details
        ptbreak.stored <- names(all.pseudotime.breakpoint.details)
    }
    while (length(tips) >= 2) {
        if (all(is.na(object@tree$segment.divergence$pseudotime.breakpoint))) {
            object@tree$segment.divergence$pseudotime.breakpoint <- 0
        }
        fuse.id <- which.max(object@tree$segment.divergence$pseudotime.breakpoint)
        seg.1 <- object@tree$segment.divergence[fuse.id, "seg.1"]
        seg.2 <- object@tree$segment.divergence[fuse.id, "seg.2"]
        pt.break <- object@tree$segment.divergence[fuse.id, "pseudotime.breakpoint"]
        if (verbose)
            print(paste("Joining segments", seg.1, "and", seg.2,
                "at pseudotime", round(pt.break, digits = 3),
                "to create segment", seg.add))
        if (is.null(object@tree$segment.joins)) {
            object@tree$segment.joins <- data.frame(child.1 = seg.1,
                child.2 = seg.2, parent = as.character(seg.add),
                pseudotime = pt.break, stringsAsFactors = F)
        }
        else {
            object@tree$segment.joins <- rbind(object@tree$segment.joins,
                c(seg.1, seg.2, as.character(seg.add), pt.break))
        }
        object@tree$segment.pseudotime.limits[seg.add, ] <- c(min(object@tree$segment.pseudotime.limits[c(seg.1,
            seg.2), "start"]), pt.break)
        object@tree$segment.pseudotime.limits[c(seg.1, seg.2),
            "start"] <- pt.break
        children.of.branchpoint <- unique(c(segChildrenAll(object,
            seg.1, include.self = T, original.joins = F, format = "binary"),
            segChildrenAll(object, seg.2, include.self = T, original.joins = F,
                format = "binary")))
        if (use.only.original.tips)
            children.of.branchpoint <- intersect(children.of.branchpoint,
                object@tree$tips)
        if (weighted.fusion) {
            weights <- tip.size[children.of.branchpoint]
            weights <- weights/sum(weights)
            object@diff.data[, paste0("visitfreq.raw.", seg.add)] <- future_apply(object@diff.data[,
                paste0("visitfreq.raw.", children.of.branchpoint)],
                1, weighted.mean, w = weights)
        }
        else {
            object@diff.data[, paste0("visitfreq.raw.", seg.add)] <- future_apply(object@diff.data[,
                paste0("visitfreq.raw.", children.of.branchpoint)],
                1, mean)
        }
        object@diff.data[, paste0("visitfreq.log.", seg.add)] <- log10(object@diff.data[,
            paste0("visitfreq.raw.", seg.add)] + 1)
        object@tree$cells.in.segments[, seg.add] <- future_apply(object@tree$cells.in.segments[,
            c(seg.1, seg.2)], 1, any)
        tips <- setdiff(c(tips, as.character(seg.add)), c(seg.1,
            seg.2))
        if (length(tips) >= 2) {
            object <- allSegmentDivergenceByPseudotime(object,
                pseudotime = pseudotime, divergence.method = divergence.method,
                segments = tips, pseudotime.cuts = cells.per.pseudotime.bin,
                window.size = bins.per.pseudotime.window, minimum.visits = minimum.visits,
                visit.threshold = visit.threshold, p.thresh = p.thresh,
                breakpoint.decision.plots = save.breakpoint.plots,
                cache = T, verbose = verbose)
            if (save.all.breakpoint.info) {
                new.breakpoints <- setdiff(names(object@tree$pseudotime.breakpoint.details),
                  ptbreak.stored)
                all.pseudotime.breakpoint.details <- unlist(list(all.pseudotime.breakpoint.details,
                  object@tree$pseudotime.breakpoint.details[new.breakpoints]),
                  recursive = F)
                ptbreak.stored <- c(ptbreak.stored, new.breakpoints)
            }
        }
        seg.add <- as.character(as.numeric(seg.add) + 1)
    }
    if (save.all.breakpoint.info)
        object@tree$pseudotime.breakpoint.details <- all.pseudotime.breakpoint.details
    else object@tree$pseudotime.breakpoint.details <- NULL
    object@tree$segments <- as.character(sort(as.numeric(unique(unlist(object@tree$segment.joins[,
        c("child.1", "child.2", "parent")])))))
    object@tree$segment.joins.initial <- object@tree$segment.joins
    if (verbose)
        print("Assigning cells to segments.")
    object <- assignCellsToSegments(object, pseudotime, verbose)
    object <- reformatSegmentJoins(object)
    object <- reformatSegmentJoins(object, segment.joins.initial = T)
    if (verbose)
        print("Collapsing short segments.")
    object <- collapseShortSegments(object, min.cells.per.segment = min.cells.per.segment,
        min.pseudotime.per.segment = min.pseudotime.per.segment)
    if (verbose)
        print("Removing singleton segments.")
    object <- removeUnitarySegments(object)
    if (verbose)
        print("Reassigning cells to segments.")
    object <- assignCellsToSegments(object, pseudotime, verbose)
    object@diff.data <- object@diff.data[unlist(object@tree$cells.in.segment),
        ]
    if (verbose)
        print("Assigning cells to nodes.")
    object <- assignCellsToNodes(object, node.size = dendro.node.size,
        pseudotime = pseudotime)
    if (verbose)
        print("Laying out tree.")
    object <- treeLayoutDendrogram(object)
    object <- treeLayoutElaborate(object)
    if (verbose)
        print("Adding cells to tree.")
    object <- treeLayoutCells(object = object, pseudotime = pseudotime,
        jitter = dendro.cell.jitter, jitter.push = dendro.cell.dist.to.tree)
    return(object)
}
environment(XY_buildTree) <- asNamespace('URD')

XY_aucprTestAlongTree <- function (object, pseudotime, tips, log.effect.size = 0.25, auc.factor = 1,
    max.auc.threshold = 1, frac.must.express = 0.1, frac.min.diff = 0.1,
    genes.use = NULL, root = NULL, segs.to.skip = NULL, only.return.global = F,
    must.beat.sibs = 0.5, report.debug = F)
{
    require(future)
    require(future.apply)
    options(future.globals.maxSize = 3000 * 1024^2)
    plan("multiprocess", workers = 8)
    plan()
    tips <- translateSegmentNames(object, tips)
    current.tip <- tips
    tip.list <- c()
    markers <- list()
    anti.markers <- list()
    cells.1.all <- c()
    cells.2.all <- c()
    stats <- data.frame(stringsAsFactors = F)
    while (length(current.tip) > 0) {
        parent.to.tip <- segParent(object, current.tip)
        if (!(current.tip %in% segs.to.skip)) {
            opposing.tips <- segSiblings(object, current.tip,
                include.self = F)
            for (opposing.tip in opposing.tips) {
                segs.to.consider <- c(opposing.tip, segChildrenAll(object,
                  opposing.tip))
                if (length(segs.to.consider) > 0) {
                  tip.list <- c(tip.list, current.tip)
                  cells.1 <- unique(unlist(object@tree$cells.in.segment[current.tip]))
                  cells.2 <- unique(unlist(object@tree$cells.in.segment[segs.to.consider]))
                  pt.limits <- object@tree$segment.pseudotime.limits[current.tip,
                    ]
                  cells.2 <- cells.2[which(object@pseudotime[cells.2,
                    pseudotime] >= as.numeric(pt.limits[1]) &
                    object@pseudotime[cells.2, pseudotime] <=
                      as.numeric(pt.limits[2]))]
                  cells.1.all <- c(cells.1.all, cells.1)
                  cells.2.all <- c(cells.2.all, cells.2)
                  if (report.debug) {
                    n.1 <- length(cells.1)
                    n.2 <- length(cells.2)
                    pt.1.mean <- mean(object@pseudotime[cells.1,
                      pseudotime])
                    pt.2.mean <- mean(object@pseudotime[cells.2,
                      pseudotime])
                    genes.1.mean <- mean(object@meta[cells.1,
                      "n.Genes"])
                    genes.2.mean <- mean(object@meta[cells.2,
                      "n.Genes"])
                    trans.1.mean <- mean(object@meta[cells.1,
                      "n.Trans"])
                    trans.2.mean <- mean(object@meta[cells.2,
                      "n.Trans"])
                    pt.1.median <- median(object@pseudotime[cells.1,
                      pseudotime])
                    pt.2.median <- median(object@pseudotime[cells.2,
                      pseudotime])
                    genes.1.median <- median(object@meta[cells.1,
                      "n.Genes"])
                    genes.2.median <- median(object@meta[cells.2,
                      "n.Genes"])
                    trans.1.median <- median(object@meta[cells.1,
                      "n.Trans"])
                    trans.2.median <- median(object@meta[cells.2,
                      "n.Trans"])
                    stats <- rbind(stats, c(n.1, n.2, pt.1.mean,
                      pt.2.mean, pt.1.median, pt.2.median, genes.1.mean,
                      genes.2.mean, genes.1.median, genes.2.median,
                      trans.1.mean, trans.2.mean, trans.1.median,
                      trans.2.median))
                  }
                  these.markers <- markersAUCPR(object = object,
                    cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use,
                    effect.size = log.effect.size, frac.must.express = frac.must.express,
                    frac.min.diff = frac.min.diff, auc.factor = auc.factor,
                    max.auc.threshold = max.auc.threshold)
                  these.anti.markers <- markersAUCPR(object = object,
                    cells.1 = cells.2, cells.2 = cells.1, genes.use = genes.use,
                    effect.size = log.effect.size, frac.must.express = frac.must.express,
                    frac.min.diff = frac.min.diff, auc.factor = auc.factor,
                    max.auc.threshold = max.auc.threshold)
                  markers[[(length(markers) + 1)]] <- these.markers[these.markers$exp.fc >
                    0, ]
                  anti.markers[[(length(anti.markers) + 1)]] <- these.anti.markers[these.anti.markers$exp.fc >
                    0, ]
                  names(markers)[length(markers)] <- paste0(current.tip,
                    "-", opposing.tip)
                  names(anti.markers)[length(anti.markers)] <- paste0(current.tip,
                    "-", opposing.tip)
                }
            }
        }
        current.tip <- parent.to.tip
        if (!is.null(root) && length(parent.to.tip) > 0 && parent.to.tip ==
            root) {
            current.tip <- NULL
        }
    }
    depleted.by.segment <- future_lapply(unique(tip.list), function(tip) {
        unique(unlist(future_lapply(anti.markers[0:(min(which(tip.list ==
            tip)) - 1)], function(markers.downstream) {
            rownames(markers.downstream)
        })))
    })
    names(depleted.by.segment) <- unique(tip.list)
    markers.2 <- future_lapply(1:length(tip.list), function(x) {
        tm <- markers[[x]]
        tm[setdiff(rownames(tm), depleted.by.segment[[tip.list[x]]]),
            ]
    })
    names(markers.2) <- names(markers)
    if (!is.null(root) && root %in% tip.list) {
        keep.markers.until <- which(tip.list == root) - 1
        markers.2 <- markers.2[1:keep.markers.until]
        tip.list <- tip.list[1:keep.markers.until]
    }
    markers.3 <- markers.2
    for (vs in unique(tip.list)) {
        sibs.use <- which(tip.list == vs)
        sib.table <- table(unlist(future_lapply(markers.2[sibs.use],
            rownames)))
        markers.keep <- names(which(sib.table >= (length(sibs.use) *
            must.beat.sibs)))
        for (this.sib in sibs.use) {
            markers.3[[this.sib]] <- markers.3[[this.sib]][intersect(rownames(markers.3[[this.sib]]),
                markers.keep), ]
        }
    }
    markers.remain <- unique(unlist(future_lapply(markers.3, rownames)))
    if (only.return.global) {
        markers.summary <- markersAUCPR(object = object, cells.1 = cells.1.all,
            cells.2 = cells.2.all, genes.use = markers.remain,
            frac.must.express = frac.must.express, effect.size = log.effect.size,
            frac.min.diff = frac.min.diff, auc.factor = auc.factor,
            max.auc.threshold = max.auc.threshold)
    }
    else {
        markers.summary <- markersAUCPR(object = object, cells.1 = cells.1.all,
            cells.2 = cells.2.all, genes.use = markers.remain,
            frac.must.express = 0, effect.size = 0, frac.min.diff = 0)
    }
    names(markers.summary) <- c("AUCPR.all", "AUCPR.ratio.all",
        "expfc.all", "posFrac_lineage", "posFrac_rest", "nTrans_lineage",
        "nTrans_rest")
    markers.max.segment <- c()
    markers.max <- future_lapply(rownames(markers.summary), function(marker) {
        id <- which.max(unlist(future_lapply(markers.3, function(x) x[marker,
            "AUCPR.ratio"])))
        to.return <- markers.3[[id]][marker, ]
        markers.max.segment <<- c(markers.max.segment, names(markers.3)[id])
        names(to.return) <- c("AUCPR_maxBranch", "AUCPR.ratio.maxBranch",
            "expfc.maxBranch", "posFrac_maxBranch", "posFrac_opposingMaxBranch",
            "nTrans_maxBranch", "nTrans_opposingMaxBranch")
        return(to.return)
    })
    markers.max <- do.call(what = "rbind", markers.max)
    markers.min.segment <- c()
    markers.min <- future_lapply(rownames(markers.summary), function(marker) {
        id <- which.min(unlist(future_lapply(markers.3, function(x) x[marker,
            "AUCPR.ratio"])))
        to.return <- markers.3[[id]][marker, ]
        markers.min.segment <<- c(markers.min.segment, names(markers.3)[id])
        names(to.return) <- c("AUCPR_minBranch", "AUCPR.ratio.minBranch",
            "expfc.minBranch", "posFrac_minBranch", "posFrac_opposingMinBranch",
            "nTrans_minBranch", "nTrans_opposingMinBranch")
        return(to.return)
    })
    markers.min <- do.call(what = "rbind", markers.min)
    markers.all <- cbind(markers.summary, markers.max, markers.min)
    markers.all$segment.maxBranch <- markers.max.segment
    markers.all$segment.minBranch <- markers.min.segment
    if (report.debug) {
        names(stats) <- c("n.1", "n.2", "pt.1.mean", "pt.2.mean",
            "pt.1.median", "pt.2.median", "genes.1.mean", "genes.2.mean",
            "genes.1.median", "genes.2.median", "trans.1.mean",
            "trans.2.mean", "trans.1.median", "trans.2.median")
        return(list(diff.exp = markers.all, stats = stats, marker.chain = markers))
    }
    else {
        return(markers.all)
    }
}
environment(XY_aucprTestAlongTree) <- asNamespace('URD')



XY_filter.genes.by.cluster.expression <- function (emat, clusters, min.max.cluster.average = 0.1){
    require(future)
    require(future.apply)
    options(future.globals.maxSize = 3000 * 1024^2)
    plan("multiprocess", workers = 8)
    plan()
    if (!any(colnames(emat) %in% names(clusters)))
        stop("provided clusters do not cover any of the emat cells!")
    vc <- intersect(colnames(emat), names(clusters))
    cl.emax <- future_apply(do.call(cbind, future_tapply(vc, as.factor(clusters[vc]),
        function(ii) Matrix::rowMeans(emat[, ii]))), 1, max)
    vi <- cl.emax > min.max.cluster.average
    emat[vi, ]
}
environment(XY_aucprTestAlongTree) <- asNamespace('velocyto.R')


XY_gene.relative.velocity.estimates <- function (emat, nmat, deltaT = 1, smat = NULL, steady.state.cells = colnames(emat),
    kCells = 10, cellKNN = NULL, kGenes = 1, old.fit = NULL,
    mult = 1000, min.nmat.smat.correlation = 0.05, min.nmat.emat.correlation = 0.05,
    min.nmat.emat.slope = 0.05, zero.offset = FALSE, deltaT2 = 1,
    fit.quantile = NULL, diagonal.quantiles = FALSE, show.gene = NULL,
    do.par = TRUE, cell.dist = NULL, emat.size = NULL, nmat.size = NULL,
    cell.emb = NULL, cell.colors = NULL, expression.gradient = NULL,
    residual.gradient = NULL, n.cores = defaultNCores(), verbose = TRUE)
{
    require(future)
    require(future.apply)
    options(future.globals.maxSize = 3000 * 1024^2)
    plan("multiprocess", workers = 8)
    plan()
    if (!all(colnames(emat) == colnames(nmat)))
        stop("emat and nmat must have the same columns (cells)")
    if (!is.null(smat)) {
        if (!all(colnames(emat) == colnames(smat)))
            stop("smat must have the same columns (cells) as emat")
    }
    resl <- list()
    vg <- intersect(rownames(emat), rownames(nmat))
    if (is.null(smat)) {
        emat <- emat[vg, ]
        nmat <- nmat[vg, ]
    }
    else {
        vg <- intersect(vg, rownames(smat))
        emat <- emat[vg, ]
        nmat <- nmat[vg, ]
        smat <- smat[vg, ]
    }
    if (!is.null(show.gene)) {
        if (!show.gene %in% rownames(emat)) {
            stop(paste("gene", show.gene, "is not present in the filtered expression matrices"))
        }
    }
    pcount <- 1
    if (!is.null(cell.dist)) {
        if (class(cell.dist) != "dist") {
            stop("cell.dist must be of a class dist")
        }
        if (!all(labels(cell.dist) == colnames(emat))) {
            cat("matching cells between cell.dist and emat/nmat ... ")
            cell.dist <- as.matrix(cell.dist)
            cn <- intersect(colnames(emat), colnames(cell.dist))
            cell.dist <- as.dist(cell.dist[cn, cn])
            emat <- emat[, cn]
            nmat <- nmat[, cn]
            if (!is.null(smat)) {
                smat <- smat[, cn]
            }
            cat("done\n")
        }
    }
    if (is.null(emat.size)) {
        emat.size <- Matrix::colSums(emat)
    }
    if (is.null(nmat.size)) {
        nmat.size <- Matrix::colSums(nmat)
    }
    emat.cs <- emat.size[colnames(emat)]/mult
    nmat.cs <- nmat.size[colnames(nmat)]/mult
    emat.log.norm <- log(as.matrix(t(t(emat)/emat.cs)) + pcount)
    if (!is.null(old.fit)) {
        cellKNN <- old.fit[["cellKNN"]]
    }
    knn.maxl <- 100
    if (kCells > 1) {
        if (is.null(cellKNN)) {
            cat("calculating cell knn ... ")
            if (is.null(cell.dist)) {
                cellKNN <- balancedKNN(emat.log.norm, kCells,
                  kCells * knn.maxl, n.threads = n.cores)
            }
            else {
                cellKNN <- balancedKNN(emat.log.norm, kCells,
                  kCells * knn.maxl, n.threads = n.cores, dist = cell.dist)
            }
            diag(cellKNN) <- 1
            resl$cellKNN <- cellKNN
            cat("done\n")
        }
        rm(emat.log.norm)
        cat("calculating convolved matrices ... ")
        conv.emat <- emat %*% cellKNN[colnames(emat), colnames(emat)]
        conv.nmat <- nmat %*% cellKNN[colnames(nmat), colnames(nmat)]
        conv.emat.cs <- (emat.cs %*% cellKNN[colnames(emat),
            colnames(emat)])[1, ]
        conv.nmat.cs <- (nmat.cs %*% cellKNN[colnames(nmat),
            colnames(nmat)])[1, ]
        cat("done\n")
    }
    else {
        conv.emat <- emat
        conv.nmat <- nmat
        cellKNN <- NULL
        conv.emat.cs <- emat.cs
        conv.nmat.cs <- nmat.cs
    }
    conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
    conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
    emat.norm <- t(t(emat)/emat.cs)
    nmat.norm <- t(t(nmat)/nmat.cs)
    if (kGenes > 1) {
        if (!is.null(old.fit) && !is.null(old.fit$geneKNN)) {
            geneKNN <- old.fit$geneKNN
        }
        else {
            cat("gene kNN ... ")
            geneKNN <- balancedKNN(t(log(as.matrix(conv.emat.norm) +
                pcount)), kGenes, kGenes * 1200, n.threads = n.cores)
            diag(geneKNN) <- 1
        }
        resl$geneKNN <- geneKNN
        cat("scaling gene weights ... ")
        gt <- rowSums(conv.emat.norm)
        scaledGeneKNN <- t(future_apply(geneKNN, 2, function(ii) pmin(1,
            median(gt[which(ii > 0)])/gt) * ii))
        cat("convolving matrices ... ")
        conv.emat.norm <- scaledGeneKNN %*% conv.emat.norm
        conv.nmat.norm <- scaledGeneKNN %*% conv.nmat.norm
        cat("done\n")
    }
    if (!is.null(smat)) {
        if (kCells > 1) {
            conv.smat <- smat %*% cellKNN[colnames(smat), colnames(smat)]
        }
        else {
            conv.smat <- smat
        }
        conv.smat.cs <- Matrix::colSums(conv.smat)/mult
        conv.smat.norm <- t(t(conv.smat)/conv.smat.cs)
        if (kGenes > 1) {
            conv.smat.norm <- scaledGeneKNN %*% conv.smat.norm
        }
        if (is.null(old.fit)) {
            cat("fitting smat-based offsets ... ")
            sfit <- data.frame(do.call(rbind, parallel::mclapply(sn(rownames(conv.emat.norm)),
                function(gn) {
                  df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                    e = (conv.emat.norm[gn, steady.state.cells]),
                    s = conv.smat.norm[gn, steady.state.cells])
                  sd <- lm(n ~ s, data = df)
                  r <- with(df[df$s > 0, ], cor(n, s, method = "spearman"),
                    3)
                  return(c(o = pmax(0, as.numeric(sd$coef[1])),
                    s = as.numeric(sd$coef[2]), r = r))
                }, mc.cores = n.cores, mc.preschedule = T)))
            cat("done\n")
        }
        else {
            sfit <- old.fit$sfit
        }
    }
    resl$conv.nmat.norm <- conv.nmat.norm
    resl$conv.emat.norm <- conv.emat.norm
    if (!is.null(show.gene)) {
        gn <- show.gene
        if (!is.null(cell.emb)) {
            cc <- intersect(rownames(cell.emb), colnames(conv.emat.norm))
            if (do.par) {
                par(mfrow = c(1, 4), mar = c(2.5, 2.5, 2.5, 0.5),
                  mgp = c(1.5, 0.65, 0), cex = 0.85)
            }
            plot(cell.emb[cc, ], pch = 21, col = ac(1, alpha = 0.2),
                bg = val2col(conv.emat.norm[gn, cc], gradientPalette = expression.gradient),
                cex = 0.8, xlab = "", ylab = "", main = paste(gn,
                  "s"), axes = F)
            box()
            plot(cell.emb[cc, ], pch = 21, col = ac(1, alpha = 0.2),
                bg = val2col(conv.nmat.norm[gn, cc], gradientPalette = expression.gradient),
                cex = 0.8, xlab = "", ylab = "", main = paste(gn,
                  "u"), axes = F)
            box()
        }
        do <- NULL
        if (!is.null(smat)) {
            df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                e = (conv.emat.norm[gn, steady.state.cells]),
                o = sfit[gn, "o"])
            if (zero.offset)
                df$o <- 0
        }
        else {
            df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                e = (conv.emat.norm[gn, steady.state.cells]))
            o <- 0
            df$o <- o
            if (!zero.offset) {
                zi <- df$e < 1/conv.emat.cs[steady.state.cells]
                if (any(zi)) {
                  o <- sum(df$n[zi])/(sum(zi) + 1)
                }
            }
            df$o <- o
        }
        d <- lm(n ~ e + offset(o) + 0, data = df, weights = df$e^4 +
            df$n^4)
        cell.col <- ac(rep(1, nrow(df)), alpha = 0.1)
        names(cell.col) <- rownames(df)
        if (!is.null(cell.colors)) {
            cc <- intersect(names(cell.colors), rownames(df))
            cell.col[cc] <- cell.colors[cc]
        }
        plot(df$e, df$n, pch = 21, bg = ac(cell.col, alpha = 0.3),
            col = ac(1, alpha = 0.1), cex = 0.8, xlab = "s",
            ylab = "u", main = paste(gn, "fit"))
        if (!is.null(do)) {
            abline(do, lty = 2, col = 8)
        }
        if (!is.null(fit.quantile)) {
            if (diagonal.quantiles) {
                emax <- quantile(df$e, p = 0.99)
                nmax <- quantile(df$n, p = 0.99)
                if (emax == 0)
                  emax <- max(max(df$e), 0.001)
                if (nmax == 0)
                  nmax <- max(max(df$n), 0.001)
                x <- df$e/emax + df$n/nmax
                eq <- quantile(x, p = c(fit.quantile, 1 - fit.quantile))
                if (!is.null(smat)) {
                  pw <- as.numeric(x >= eq[2])
                }
                else {
                  pw <- as.numeric(x >= eq[2] | x <= eq[1])
                }
            }
            else {
                eq <- quantile(df$e, p = c(fit.quantile, 1 -
                  fit.quantile))
                if (!is.null(smat) || zero.offset) {
                  pw <- as.numeric(df$e >= eq[2])
                }
                else {
                  pw <- as.numeric(df$e >= eq[2] | df$e <= eq[1])
                }
            }
            if (!is.null(smat) || zero.offset) {
                d <- lm(n ~ e + offset(o) + 0, data = df, weights = pw)
            }
            else {
                d <- lm(n ~ e, data = df, weights = pw)
            }
        }
        df <- df[order(df$e, decreasing = T), ]
        lines(df$e, predict(d, newdata = df), lty = 2, col = 2)
        if (!is.null(cell.emb)) {
            plot(cell.emb[cc, ], pch = 21, col = ac(1, alpha = 0.2),
                bg = val2col(resid(d)[cc], gradientPalette = residual.gradient),
                cex = 0.8, xlab = "", ylab = "", main = paste(gn,
                  "resid"), axes = F)
            box()
        }
        if (kGenes > 1) {
            return(invisible(geneKNN))
        }
        else {
            return(1)
        }
    }
    cat("fitting gamma coefficients ... ")
    if (is.null(old.fit)) {
        ko <- data.frame(do.call(rbind, parallel::mclapply(sn(rownames(conv.emat.norm)),
            function(gn) {
                if (!is.null(smat)) {
                  df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                    e = (conv.emat.norm[gn, steady.state.cells]),
                    o = sfit[gn, "o"])
                  if (zero.offset)
                    df$o <- 0
                }
                else {
                  df <- data.frame(n = (conv.nmat.norm[gn, steady.state.cells]),
                    e = (conv.emat.norm[gn, steady.state.cells]))
                  o <- 0
                  if (!zero.offset) {
                    zi <- df$e < 1/conv.emat.cs[steady.state.cells]
                    if (any(zi)) {
                      o <- sum(df$n[zi])/(sum(zi) + 1)
                    }
                  }
                  df$o <- o
                }
                if (is.null(fit.quantile)) {
                  d <- lm(n ~ e + offset(o) + 0, data = df, weights = df$e^4 +
                    df$n^4)
                  return(c(o = df$o[1], g = as.numeric(coef(d)[1]),
                    r = cor(df$e, df$n, method = "spearman")))
                }
                else {
                  if (diagonal.quantiles) {
                    emax <- quantile(df$e, p = 0.99)
                    nmax <- quantile(df$n, p = 0.99)
                    if (emax == 0)
                      emax <- max(max(df$e), 0.001)
                    if (nmax == 0)
                      nmax <- max(max(df$n), 0.001)
                    x <- df$e/emax + df$n/nmax
                    eq <- quantile(x, p = c(fit.quantile, 1 -
                      fit.quantile))
                    if (!is.null(smat)) {
                      pw <- as.numeric(x >= eq[2])
                    }
                    else {
                      pw <- as.numeric(x >= eq[2] | x <= eq[1])
                    }
                  }
                  else {
                    eq <- quantile(df$e, p = c(fit.quantile,
                      1 - fit.quantile))
                    if (!is.null(smat) || zero.offset) {
                      pw <- as.numeric(df$e >= eq[2])
                    }
                    else {
                      pw <- as.numeric(df$e >= eq[2] | df$e <=
                        eq[1])
                    }
                  }
                  if (!is.null(smat) || zero.offset) {
                    d <- lm(n ~ e + offset(o) + 0, data = df,
                      weights = pw)
                    return(c(o = df$o[1], g = as.numeric(coef(d)[1]),
                      r = cor(df$e, df$n, method = "spearman")))
                  }
                  else {
                    d <- lm(n ~ e, data = df, weights = pw)
                    return(c(o = as.numeric(coef(d)[1]), g = as.numeric(coef(d)[2]),
                      r = cor(df$e, df$n, method = "spearman")))
                  }
                }
            }, mc.cores = n.cores, mc.preschedule = T)))
        ko <- na.omit(ko)
        cat("done. succesfful fit for", nrow(ko), "genes\n")
    }
    else {
        full.ko <- ko <- na.omit(old.fit$ko)
    }
    if (!is.null(smat)) {
        sfit <- na.omit(sfit)
        ko <- ko[rownames(ko) %in% rownames(sfit), ]
        vi <- sfit$r > min.nmat.smat.correlation
        ko <- ko[vi, ]
        if (!all(vi))
            cat("filtered out", sum(!vi), "out of", length(vi),
                "genes due to low nmat-smat correlation\n")
    }
    full.ko <- ko
    vi <- ko$r > min.nmat.emat.correlation
    if (!all(vi))
        cat("filtered out", sum(!vi), "out of", length(vi), "genes due to low nmat-emat correlation\n")
    ko <- ko[vi, ]
    vi <- ko$g > min.nmat.emat.slope
    if (!all(vi))
        cat("filtered out", sum(!vi), "out of", length(vi), "genes due to low nmat-emat slope\n")
    ko <- ko[vi, ]
    gamma <- ko$g
    offset <- ko$o
    names(gamma) <- names(offset) <- rownames(ko)
    cat("calculating RNA velocity shift ... ")
    if (kGenes > 1) {
        npred <- gamma * conv.emat.norm[names(gamma), ] + ko$o
        npred[npred < 0] <- 0
        mval <- log2(conv.nmat.norm[names(gamma), ] + pcount) -
            log2(npred + pcount)
        resl$mval <- mval
        conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
        conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
        cat("re-estimating gamma of individual genes ... ")
        am <- conv.nmat.norm[rownames(mval), ] - offset
        am[am < 0] <- 0
        fm <- log2(am) - mval - log2(conv.emat.norm[rownames(mval),
            ])
        wm <- is.finite(fm)
        fm[!is.finite(fm)] <- 0
        gammaA <- 2^(rowSums(fm * wm)/rowSums(wm))
        gammaA <- gammaA[is.finite(gammaA)]
        gamma <- gammaA
        cat("done\n")
        cat("calculating RNA velocity shift ... ")
        deltaE <- t.get.projected.delta.from.log2ratio(em = conv.emat.norm,
            gamma = gamma, r = mval, delta = deltaT)
    }
    else {
        deltaE <- t.get.projected.delta(conv.emat.norm, conv.nmat.norm,
            gamma, offset = offset, delta = deltaT)
    }
    resl$gamma <- gamma
    cat("done\n")
    cat("calculating extrapolated cell state ... ")
    emat.norm <- emat[rownames(emat) %in% rownames(deltaE), ]
    emat.sz <- emat.cs
    emat.norm <- t(t(emat.norm)/(emat.sz))
    emn <- t.get.projected.cell2(emat.norm, emat.sz, as.matrix(deltaE),
        mult = mult, delta = deltaT2)
    cat("done\n")
    full.ko$valid <- rownames(full.ko) %in% rownames(ko)
    resl <- c(resl, list(projected = emn, current = emat.norm,
        deltaE = deltaE, deltaT = deltaT, ko = full.ko, mult = mult,
        kCells = kCells))
    if (!is.null(smat)) {
        resl$sfit <- sfit
    }
    return(resl)
}
environment(XY_gene.relative.velocity.estimates) <- asNamespace('velocyto.R')

XY_show.velocity.on.embedding.cor <- function (emb, vel, n = 100, cell.colors = NULL, corr.sigma = 0.05,
    show.grid.flow = FALSE, grid.n = 20, grid.sd = NULL, min.grid.cell.mass = 1,
    min.arrow.size = NULL, arrow.scale = 1, max.grid.arrow.length = NULL,
    fixed.arrow.length = FALSE, plot.grid.points = FALSE, scale = "log",
    nPcs = NA, arrow.lwd = 1, xlab = "", ylab = "", n.cores = defaultNCores(),
    do.par = T, show.cell = NULL, cell.border.alpha = 0.3, cc = NULL,
    return.details = FALSE, expression.scaling = FALSE, ...)
{
    randomize <- FALSE
    if (do.par)
        par(mfrow = c(1, 1), mar = c(3.5, 3.5, 2.5, 1.5), mgp = c(2,
            0.65, 0), cex = 0.85)
    celcol <- "white"
    if (is.null(show.cell)) {
        celcol <- cell.colors[rownames(emb)]
    }
    plot(emb, bg = celcol, pch = 21, col = ac(1, alpha = cell.border.alpha),
        xlab = xlab, ylab = ylab, ...)
    em <- as.matrix(vel$current)
    ccells <- intersect(rownames(emb), colnames(em))
    em <- em[, ccells]
    emb <- emb[ccells, ]
    nd <- as.matrix(vel$deltaE[, ccells])
    cgenes <- intersect(rownames(em), rownames(nd))
    nd <- nd[cgenes, ]
    em <- em[cgenes, ]
    if (randomize) {
        nd <- t(future_apply(nd, 1, function(x) (rbinom(length(x), 1,
            0.5) * 2 - 1) * abs(sample(x))))
    }
    if (is.null(cc)) {
        cat("delta projections ... ")
        if (scale == "log") {
            cat("log ")
            cc <- colDeltaCorLog10(em, (log10(abs(nd) + 1) *
                sign(nd)), nthreads = n.cores)
        }
        else if (scale == "sqrt") {
            cat("sqrt ")
            cc <- colDeltaCorSqrt(em, (sqrt(abs(nd)) * sign(nd)),
                nthreads = n.cores)
        }
        else if (scale == "rank") {
            cat("rank ")
            cc <- colDeltaCor((future_apply(em, 2, rank)), (future_apply(nd,
                2, rank)), nthreads = n.cores)
        }
        else {
            cat("linear ")
            cc <- colDeltaCor(em, nd, nthreads = n.cores)
        }
        colnames(cc) <- rownames(cc) <- colnames(em)
        diag(cc) <- 0
    }
    cat("knn ... ")
    if (n > nrow(cc)) {
        n <- nrow(cc)
    }
    emb.knn <- balancedKNN(t(emb), k = n, maxl = nrow(emb), dist = "euclidean",
        n.threads = n.cores)
    diag(emb.knn) <- 1
    cat("transition probs ... ")
    tp <- exp(cc/corr.sigma) * emb.knn
    tp <- t(t(tp)/Matrix::colSums(tp))
    tp <- as(tp, "dgCMatrix")
    cat("done\n")
    if (!is.null(show.cell)) {
        i <- match(show.cell, rownames(emb))
        if (is.na(i))
            stop(paste("specified cell", i, "is not in the embedding"))
        points(emb, pch = 19, col = ac(val2col(tp[rownames(emb),
            show.cell], gradient.range.quantile = 1), alpha = 0.5))
        points(emb[show.cell, 1], emb[show.cell, 2], pch = 3,
            cex = 1, col = 1)
        di <- t(t(emb) - emb[i, ])
        di <- di/sqrt(Matrix::rowSums(di^2)) * arrow.scale
        di[i, ] <- 0
        dir <- Matrix::colSums(di * tp[, i])
        dic <- Matrix::colSums(di * (tp[, i] > 0)/sum(tp[, i] >
            0))
        dia <- dir - dic
        suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i],
            2], emb[colnames(em)[i], 1] + dic[1], emb[colnames(em)[i],
            2] + dic[2], length = 0.05, lwd = 1, col = "blue"))
        suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i],
            2], emb[colnames(em)[i], 1] + dir[1], emb[colnames(em)[i],
            2] + dir[2], length = 0.05, lwd = 1, col = "red"))
        suppressWarnings(arrows(emb[colnames(em)[i], 1] + dic[1],
            emb[colnames(em)[i], 2] + dic[2], emb[colnames(em)[i],
                1] + dir[1], emb[colnames(em)[i], 2] + dir[2],
            length = 0.05, lwd = 1, lty = 1, col = "grey50"))
        suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i],
            2], emb[colnames(em)[i], 1] + dia[1], emb[colnames(em)[i],
            2] + dia[2], length = 0.05, lwd = 1, col = "black"))
    }
    else {
        cat("calculating arrows ... ")
        arsd <- data.frame(t(embArrows(emb, tp, arrow.scale,
            n.cores)))
        rownames(arsd) <- rownames(emb)
        if (expression.scaling) {
            tpb <- tp > 0
            tpb <- t(t(tpb)/colSums(tpb))
            es <- as.matrix(em %*% tp) - as.matrix(em %*% as.matrix(tpb))
            pl <- pmin(1, pmax(0, future_apply(as.matrix(vel$deltaE[,
                colnames(es)]) * es, 2, sum)/sqrt(colSums(es *
                es))))
            arsd <- arsd * pl
        }
        ars <- data.frame(cbind(emb, emb + arsd))
        colnames(ars) <- c("x0", "y0", "x1", "y1")
        colnames(arsd) <- c("xd", "yd")
        rownames(ars) <- rownames(emb)
        cat("done\n")
        if (show.grid.flow) {
            cat("grid estimates ... ")
            rx <- range(c(range(ars$x0), range(ars$x1)))
            ry <- range(c(range(ars$y0), range(ars$y1)))
            gx <- seq(rx[1], rx[2], length.out = grid.n)
            gy <- seq(ry[1], ry[2], length.out = grid.n)
            if (is.null(grid.sd)) {
                grid.sd <- sqrt((gx[2] - gx[1])^2 + (gy[2] -
                  gy[1])^2)/2
                cat("grid.sd=", grid.sd, " ")
            }
            if (is.null(min.arrow.size)) {
                min.arrow.size <- sqrt((gx[2] - gx[1])^2 + (gy[2] -
                  gy[1])^2) * 0.01
                cat("min.arrow.size=", min.arrow.size, " ")
            }
            if (is.null(max.grid.arrow.length)) {
                max.grid.arrow.length <- sqrt(sum((par("pin")/c(length(gx),
                  length(gy)))^2)) * 0.25
                cat("max.grid.arrow.length=", max.grid.arrow.length,
                  " ")
            }
            garrows <- do.call(rbind, future_lapply(gx, function(x) {
                cd <- sqrt(outer(emb[, 2], -gy, "+")^2 + (x -
                  emb[, 1])^2)
                cw <- dnorm(cd, sd = grid.sd)
                gw <- Matrix::colSums(cw)
                cws <- pmax(1, Matrix::colSums(cw))
                gxd <- Matrix::colSums(cw * arsd$xd)/cws
                gyd <- Matrix::colSums(cw * arsd$yd)/cws
                al <- sqrt(gxd^2 + gyd^2)
                vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
                cbind(rep(x, sum(vg)), gy[vg], x + gxd[vg], gy[vg] +
                  gyd[vg])
            }))
            colnames(garrows) <- c("x0", "y0", "x1", "y1")
            if (fixed.arrow.length) {
                suppressWarnings(arrows(garrows[, 1], garrows[,
                  2], garrows[, 3], garrows[, 4], length = 0.05,
                  lwd = arrow.lwd))
            }
            else {
                alen <- pmin(max.grid.arrow.length, sqrt(((garrows[,
                  3] - garrows[, 1]) * par("pin")[1]/diff(par("usr")[c(1,
                  2)]))^2 + ((garrows[, 4] - garrows[, 2]) *
                  par("pin")[2]/diff(par("usr")[c(3, 4)]))^2))
                suppressWarnings(future_lapply(1:nrow(garrows), function(i) arrows(garrows[i,
                  1], garrows[i, 2], garrows[i, 3], garrows[i,
                  4], length = alen[i], lwd = arrow.lwd)))
            }
            if (plot.grid.points)
                points(rep(gx, each = length(gy)), rep(gy, length(gx)),
                  pch = ".", cex = 0.1, col = ac(1, alpha = 0.4))
            cat("done\n")
            if (return.details) {
                cat("expression shifts .")
                scale.int <- switch(scale, log = 2, sqrt = 3,
                  1)
                if (!expression.scaling) {
                  tpb <- tp > 0
                  tpb <- t(t(tpb)/colSums(tpb))
                  es <- as.matrix(em %*% tp) - as.matrix(em %*%
                    as.matrix(tpb))
                }
                cat(".")
                gs <- do.call(cbind, parallel::mclapply(gx, function(x) {
                  cd <- sqrt(outer(emb[, 2], -gy, "+")^2 + (x -
                    emb[, 1])^2)
                  cw <- dnorm(cd, sd = grid.sd)
                  gw <- Matrix::colSums(cw)
                  cws <- pmax(1, Matrix::colSums(cw))
                  cw <- t(t(cw)/cws)
                  gxd <- Matrix::colSums(cw * arsd$xd)
                  gyd <- Matrix::colSums(cw * arsd$yd)
                  al <- sqrt(gxd^2 + gyd^2)
                  vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
                  if (any(vg)) {
                    z <- es %*% cw[, vg]
                  }
                  else {
                    NULL
                  }
                }, mc.cores = n.cores, mc.preschedule = T))
                if (scale == "log") {
                  nd <- (log10(abs(nd) + 1) * sign(nd))
                }
                else if (scale == "sqrt") {
                  nd <- (sqrt(abs(nd)) * sign(nd))
                }
                cat(".")
                gv <- do.call(cbind, parallel::mclapply(gx, function(x) {
                  cd <- sqrt(outer(emb[, 2], -gy, "+")^2 + (x -
                    emb[, 1])^2)
                  cw <- dnorm(cd, sd = grid.sd)
                  gw <- Matrix::colSums(cw)
                  cws <- pmax(1, Matrix::colSums(cw))
                  cw <- t(t(cw)/cws)
                  gxd <- Matrix::colSums(cw * arsd$xd)
                  gyd <- Matrix::colSums(cw * arsd$yd)
                  al <- sqrt(gxd^2 + gyd^2)
                  vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
                  if (any(vg)) {
                    z <- nd %*% cw[, vg]
                  }
                  else {
                    NULL
                  }
                }, mc.cores = n.cores, mc.preschedule = T))
                cat(". done\n")
                return(invisible(list(tp = tp, cc = cc, garrows = garrows,
                  arrows = as.matrix(ars), vel = nd, eshifts = es,
                  gvel = gv, geshifts = gs, scale = scale)))
            }
        }
        else {
            future_apply(ars, 1, function(x) {
                if (fixed.arrow.length) {
                  suppressWarnings(arrows(x[1], x[2], x[3], x[4],
                    length = 0.05, lwd = arrow.lwd))
                }
                else {
                  ali <- sqrt(((x[3] - x[1]) * par("pin")[1]/diff(par("usr")[c(1,
                    2)]))^2 + ((x[4] - x[2]) * par("pin")[2]/diff(par("usr")[c(3,
                    4)]))^2)
                  suppressWarnings(arrows(x[1], x[2], x[3], x[4],
                    length = min(0.05, ali), lwd = arrow.lwd))
                }
            })
        }
    }
    return(invisible(list(tp = tp, cc = cc)))
}
environment(XY_show.velocity.on.embedding.cor) <- asNamespace('velocyto.R')




XY_wbm_order_new_DimPlot <- function (object, dims = c(1, 2), cells = NULL, cols = NULL,
    pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL,
    shape.by = NULL, order = NULL, label = FALSE, label.size = 4,
    repel = FALSE, cells.highlight = NULL, cols.highlight = "red",
    sizes.highlight = 1, na.value = "grey50", combine = TRUE,
    ...)
{
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    reduction <- reduction %||% {
        default.reductions <- c("umap", "tsne", "pca")
        object.reductions <- FilterObjects(object = object, classes.keep = "DimReduc")
        reduc.use <- min(which(x = default.reductions %in% object.reductions))
        default.reductions[reduc.use]
    }
    cells <- cells %||% colnames(x = object)
    data <- Embeddings(object = object[[reduction]])[cells, dims]
    data <- as.data.frame(x = data)
    dims <- paste0(Key(object = object[[reduction]]), dims)
    object <- suppressMessages(expr = StashIdent(object = object,
        save.name = "ident"))
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
        if (length(x = group.by) > 1) {
            nrow <- ncol <- 1
        }
        else {
            nrow <- NULL
            ncol <- list(...)$ncol
        }
        data[, split.by] <- object[[split.by, drop = TRUE]]
    }
    data
    data <- data[order((data$new),decreasing=F),]
    plots <- lapply(X = group.by, FUN = function(x) {
        return(SingleDimPlot(data = data[, c(dims, x, split.by)],
            dims = dims, col.by = x, cols = cols, pt.size = pt.size,
            shape.by = shape.by, label = label,
            repel = repel, label.size = label.size, cells.highlight = cells.highlight,
            cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
            na.value = na.value))
    })
    if (!is.null(x = split.by)) {
        plots <- lapply(X = plots, FUN = function(x) {
            x + facet_wrap(facets = split.by, nrow = 1)
        })
    }
    if (combine) {
        plots <- CombinePlots(plots = plots, ncol = 4, ...)
    }
    return(plots)
}
environment(XY_wbm_order_new_DimPlot) <- asNamespace('Seurat')


Pesudo_FeaturePlot <- function (object, features, dims = c(1, 2), cells = NULL, cols = c("lightgrey",
      "blue"), pt.size = NULL, min.cutoff = NA, max.cutoff = NA,
      reduction = NULL, split.by = NULL, shape.by = NULL, blend = FALSE,
      blend.threshold = 0.5, order = NULL, label = FALSE, label.size = 4,
      ncol = NULL, combine = TRUE, coord.fixed = FALSE, ...)
 {
      no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",
              size = 14, margin = margin(r = 7)))
      if (is.null(reduction)) {
          default_order <- c("umap", "tsne", "pca")
          reducs <- which(default_order %in% names(object@reductions))
          reduction <- default_order[reducs[1]]
      }
      if (length(x = dims) != 2 || !is.numeric(x = dims)) {
          stop("'dims' must be a two-length integer vector")
      }
      if (blend && length(x = features) != 2) {
          stop("Blending feature plots only works with two features")
      }
      dims <- paste0(Key(object = object[[reduction]]), dims)
      cells <- cells %||% colnames(x = object)
      data <- FetchData(object = object, vars = c(dims, features),
          cells = cells)
      features <- colnames(x = data)[3:ncol(x = data)]
      min.cutoff <- mapply(FUN = function(cutoff, feature) {
          return(ifelse(test = is.na(x = cutoff), yes = min(data[,
              feature]), no = cutoff))
      }, cutoff = min.cutoff, feature = features)
      max.cutoff <- mapply(FUN = function(cutoff, feature) {
          return(ifelse(test = is.na(x = cutoff), yes = max(data[,
              feature]), no = cutoff))
      }, cutoff = max.cutoff, feature = features)
      check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
          max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
      if (length(x = check.lengths) != 1) {
          stop("There must be the same number of minimum and maximum cuttoffs as there are features")
      }
      brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
          ]$maxcolors, no = length(x = cols))
      data[, 3:ncol(x = data)] <- sapply(X = 3:ncol(x = data),
          FUN = function(index) {
              data.feature <- as.vector(x = data[, index])
              min.use <- SetQuantile(cutoff = min.cutoff[index -
                  2], data.feature)
              max.use <- SetQuantile(cutoff = max.cutoff[index -
                  2], data.feature)
              data.feature[data.feature < min.use] <- min.use
              data.feature[data.feature > max.use] <- max.use
              if (brewer.gran == 2) {
                  return(data.feature)
              }
              data.cut <- if (all(data.feature == 0)) {
                  0
              }
              else {
                  as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                    breaks = brewer.gran)))
              }
              return(data.cut)
          })
      colnames(x = data)[3:ncol(x = data)] <- features
      rownames(x = data) <- cells
      data$split <- if (is.null(x = split.by)) {
          RandomName()
      }
      else {
          switch(EXPR = split.by, ident = Idents(object = object)[cells],
              object[[split.by, drop = TRUE]][cells])
      }
      if (!is.factor(x = data$split)) {
          data$split <- factor(x = data$split)
      }
      plots <- vector(mode = "list", length = ifelse(test = blend,
          yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
      xlims <- c(min(data[, dims[1]]), max(data[,dims[1]]))
      ylims <- c(min(data[, dims[2]]), max(data[,dims[2]]))
      if (blend) {
          ncol <- 4
          color.matrix <- BlendMatrix(col.threshold = blend.threshold)
          colors <- list(color.matrix[, 1], color.matrix[1, ],
              as.vector(x = color.matrix))
      }
      for (i in 1:length(x = levels(x = data$split))) {
          ident <- levels(x = data$split)[i]
          data.plot <- data[as.character(x = data$split) == ident,
              , drop = FALSE]
          if (blend) {
              data.plot <- cbind(data.plot[, dims], BlendExpression(data = data.plot[,
                  features[1:2]]))
              features <- colnames(x = data.plot)[3:ncol(x = data.plot)]
          }
          for (j in 1:length(x = features)) {
              feature <- features[j]
              if (blend) {
                  cols.use <- as.numeric(x = as.character(x = data.plot[,
                    feature])) + 1
                  cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
              }
              else {
                  cols.use <- NULL
              }
              data.plot <- data.plot[order(data.plot[,j+2],decreasing=F,na.last=F),]
              plot <- SingleDimPlot(data = data.plot[, c(dims,
                  feature)], dims = dims, col.by = feature, pt.size = pt.size,
                  cols = cols.use, label.size = label.size) +
                  scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
                  theme_cowplot()
              if (length(x = levels(x = data$split)) > 1) {
                  plot <- plot + theme(panel.border = element_rect(fill = NA,
                    colour = "black"))
                  plot <- plot + if (i == 1) {
                    labs(title = feature)
                  }
                  else {
                    labs(title = NULL)
                  }
                  if (j == length(x = features) && !blend) {
                    suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident)) +
                      no.right)
                  }
                  if (j != 1) {
                    plot <- plot + theme(axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                      axis.title.y.left = element_blank())
                  }
                  if (i != length(x = levels(x = data$split))) {
                    plot <- plot + theme(axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                      axis.title.x = element_blank())
                  }
              }
              else {
                  plot <- plot + labs(title = feature)
              }
              if (!blend) {
                  plot <- plot + guides(color = NULL)
                  if (length(x = cols) == 1) {
                    plot <- plot + scale_color_brewer(palette = cols)
                  }
                  else if (length(x = cols) > 1) {
                    plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols,
                      guide = "colorbar"))
                  }
              }
              if (coord.fixed) {
                  plot <- plot + coord_fixed()
              }
              plot <- plot
              plots[[(length(x = features) * (i - 1)) + j]] <- plot
          }
      }
      if (blend) {
          blend.legend <- BlendMap(color.matrix = color.matrix)
          for (i in 1:length(x = levels(x = data$split))) {
              suppressMessages(expr = plots <- append(x = plots,
                  values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                    1, yes = levels(x = data$split)[i], no = "")),
                    expand = c(0, 0)) + labs(x = features[1], y = features[2],
                    title = if (i == 1) {
                      paste("Color threshold:", blend.threshold)
                    } else {
                      NULL
                    }) + no.right), after = 4 * i - 1))
          }
      }
      plots <- Filter(f = Negate(f = is.null), x = plots)
      if (combine) {
          if (is.null(x = ncol)) {
              ncol <- 2
              if (length(x = features) == 1) {
                  ncol <- 1
              }
              if (length(x = features) > 6) {
                  ncol <- 3
              }
              if (length(x = features) > 9) {
                  ncol <- 4
              }
          }
          ncol <- ifelse(test = is.null(x = split.by) || blend,
              yes = ncol, no = length(x = features))
          legend <- if (blend) {
              "none"
          }
          else {
              split.by %iff% "none"
          }
          plots <- CombinePlots(plots = plots, ncol = ncol, legend = legend,
              nrow = split.by %iff% length(x = levels(x = data$split)))
      }
      return(plots)
}
environment(Pesudo_FeaturePlot) <- asNamespace('Seurat')


XY_Pesudo_heatmap <- function (scale_data=scale_data, group=group,anno_data=anno_data,color=color,
  min_and_max_cut=min_and_max_cut,num_cluster=num_cluster,hclust_method=hclust_method,cluster_rows=TRUE){
  message("Processed data begain")
  sel_cutoff <- min(abs(range(scale_data)))
  filter_scale_data <- scale_data
  if (is.null(min_and_max_cut) | sel_cutoff < min_and_max_cut){
    filter_scale_data[filter_scale_data > sel_cutoff] <- sel_cutoff
    filter_scale_data[filter_scale_data < -sel_cutoff] <- -sel_cutoff
    } else {
      filter_scale_data[filter_scale_data > min_and_max_cut] <- min_and_max_cut
      filter_scale_data[filter_scale_data < -min_and_max_cut] <- -min_and_max_cut
    }
    message("datasets range in ",range(filter_scale_data)[1]," to ",range(filter_scale_data)[2])
    if (length(table(rownames(anno_data) %in% colnames(filter_scale_data)))==2){
      stop("No matched anno_data and filter_scale_data")
    }
    anno_data <- anno_data[colnames(filter_scale_data),]
  col_annotation = data.frame(new_anno=anno_data[group],row.names=rownames(anno_data))
  colnames(col_annotation) <- "NewAnno"
  require(pheatmap)
  if (is.null(num_cluster)){
    row_cluster <- NULL
  } else {
    row_cluster <- num_cluster
  }
  message("cluster begain")
  row_dist <- as.dist((1 - cor(Matrix::t(filter_scale_data)))/2)
  row_dist[is.na(row_dist)] <- 1
  ph <- pheatmap(filter_scale_data, useRaster = T, cluster_cols = FALSE,
    cluster_rows = cluster_rows, show_rownames = F, show_colnames = F,
    clustering_distance_rows = row_dist, clustering_method = hclust_method,
    cutree_rows = row_cluster, silent = TRUE, filename = NA,
    border_color = NA, color = rev(color))
  annotation_row <- data.frame(ColClu = factor(cutree(ph$tree_row,row_cluster)))
  annotation_row$ColClu <- paste("Clu",annotation_row$ColClu,sep="_")
  clu_sel_cols <- brewer.pal(12, "Set3")
  clu_col <- clu_sel_cols[1:row_cluster]
  names(clu_col) <- as.character(unique(annotation_row$ColClu))
  anno_sel_cols <- brewer.pal(8, "Pastel1")
  anno_col <- anno_sel_cols[1:length(unique(col_annotation$NewAnno))]
  names(anno_col) <- as.character(unique(col_annotation$NewAnno))
  ann_colors = list(ColClu = clu_col,NewAnno = anno_col)
  message("pheatmap printing start")
  ph_res <- pheatmap(filter_scale_data, useRaster = T, cluster_cols = FALSE, annotation_colors = ann_colors,
    cluster_rows = cluster_rows, show_rownames = F,
    show_colnames = F, clustering_distance_rows = row_dist,
    clustering_method = hclust_method, cutree_rows = row_cluster,
    annotation_row = annotation_row, annotation_col = col_annotation,
    treeheight_row = 20, fontsize = 6, color = rev(color),
    border_color = NA, silent = TRUE, filename = NA)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  return(ph_res)
}





XY_cnetplot.enrichResult <- function(x,showCategory = 5,foldChange   = NULL,layout = "kk",colorEdge = FALSE,circular = FALSE,node_label = "all",...) {
  node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
  if (circular) {
          layout <- "linear"
          geom_edge <- geom_edge_arc
  } else {
      geom_edge <- geom_edge_link
  }
  geneSets <- extract_geneSets(x, showCategory)
  g <- list2graph(geneSets)
  foldChange <- fc_readable(x, foldChange)
  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  if (colorEdge) {
          E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
          edge_layer <- geom_edge(aes_(color = ~category), alpha=.8)
  } else {
      edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
  }
  if (!is.null(foldChange)) {
          fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
          V(g)$color <- NA
          V(g)$color[(n+1):length(V(g))] <- fc
          palette <- fc_palette(fc)
          p <- ggraph(g, layout=layout, circular = circular) +
              edge_layer +
              geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
              scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494")
  } else {
          V(g)$color <- "#B3B3B3"
          V(g)$color[1:n] <- "#E5C494"
          p <- ggraph(g, layout=layout, circular=circular) +
              edge_layer +
              geom_node_point(aes_(color=~I(color), size=~size))
  }
  p <- p + scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
      theme_void()
  if (node_label == "category") {
      p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,])
  } else if (node_label == "gene") {
      p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),], repel=TRUE)
  } else if (node_label == "all") {
      p <- p + geom_node_text(aes_(label=~name), repel=TRUE)
  } 
return(p)
}
environment(XY_cnetplot.enrichResult) <- asNamespace('enrichplot')







XY_ImportSeuratObject <- function (seuratobj, clusters, timepoints, timepoint_order, cluster_labels){
    if (class(seuratobj)[1] == "Seurat"){
        requireNamespace("Seurat")
        } else {
            stop("Not a Seurat object. Tempora only supports importing Seurat objects at the moment. See ?Tempora::CreateTemporaObject to manually create a Tempora object from an expression matrix")
        }
        data <- GetAssayData(seuratobj,slot="data")
        cat("Extracting data...")
        metadata <- seuratobj@meta.data
        cat("\nExtracting metadata...")
        colnames(metadata)[which(colnames(metadata) == clusters)] <- "Clusters"
        colnames(metadata)[which(colnames(metadata) == timepoints)] <- "Timepoints"
        cat("\nCreating Tempora object...")
        tempora_obj <- CreateTemporaObject(as.matrix(data), meta.data = metadata,
            timepoint_order = timepoint_order, cluster_labels = cluster_labels)
        validObject(tempora_obj)
        return(tempora_obj)
    }




XY_read.loom.matrices <- function (file, engine = "hdf5r"){
  require(future.apply)
    if (engine == "h5") {
        cat("reading loom file via h5...\n")
        f <- h5::h5file(file, mode = "r")
        cells <- f["col_attrs/CellID"][]
        genes <- f["row_attrs/Gene"][]
        dl <- c(spliced = "/layers/spliced", unspliced = "/layers/unspliced",
            ambiguous = "/layers/ambiguous")
        if ("/layers/spanning" %in% h5::list.datasets(f)) {
            dl <- c(dl, c(spanning = "/layers/spanning"))
        }
        dlist <- future_lapply(dl, function(path) {
            m <- as(f[path][], "dgCMatrix")
            rownames(m) <- genes
            colnames(m) <- cells
            return(m)
        })
        h5::h5close(f)
        return(dlist)
    }
    else if (engine == "hdf5r") {
        cat("reading loom file via hdf5r...\n")
        f <- hdf5r::H5File$new(file, mode = "r")
        cells <- f[["col_attrs/CellID"]][]
        genes <- f[["row_attrs/Gene"]][]
        dl <- c(spliced = "layers/spliced", unspliced = "layers/unspliced",
            ambiguous = "layers/ambiguous")
        if ("layers/spanning" %in% hdf5r::list.datasets(f)) {
            dl <- c(dl, c(spanning = "layers/spanning"))
        }
        dlist <- future_lapply(dl, function(path) {
            m <- as(t(f[[path]][, ]), "dgCMatrix")
            rownames(m) <- genes
            colnames(m) <- cells
            return(m)
        })
        f$close_all()
        return(dlist)
    }
    else {
        warning("Unknown engine. Use hdf5r or h5 to import loom file.")
        return(list())
    }
}
environment(XY_read.loom.matrices) <- asNamespace('velocyto.R')



update_n <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
        return(showCategory)
    }

    ## geneSets <- geneInCategory(x) ## use core gene for gsea result
    n <- showCategory
    if (nrow(x) < n) {
        n <- nrow(x)
    }

    return(n)
}

extract_geneSets <- function(x, n) {
    n <- update_n(x, n)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description
    if (is.numeric(n)) {
        return(geneSets[1:n])
    }
    return(geneSets[n]) ## if n is a vector of Description
}


list2df <- function(inputList) {
    ldf <- lapply(1:length(inputList), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })

    do.call('rbind', ldf)
}
list2graph <- function(inputList) {
    x <- list2df(inputList)
    g <- graph.data.frame(x, directed=FALSE)
    return(g)
}



XY_PlotCoverage <- function (genome_gr, geneSymbol = "", wig_data = NULL, bamfiles = NULL,
    peaks.annot = NULL, label.transcripts = FALSE, wig_same_strand = TRUE,
    genome = NULL, pdf_output = FALSE, wig_data.tracknames = NULL,
    bamfile.tracknames = NULL, output_file_name = "", zoom_3UTR = FALSE,
    annotation.fontsize = NULL, axis.fontsize = NULL)
{
    GenomeInfoDb::seqlevelsStyle(genome_gr) <- "UCSC"
    idx <- which(genome_gr$gene_name == geneSymbol)
    if (length(idx) == 0) {
        warning("Could not find gene name. Please check spelling (and case)")
        return(NULL)
    }
    genome_gr <- genome_gr[idx]
    start <- min(IRanges::start(IRanges::ranges(genome_gr)))
    end <- max(IRanges::end(IRanges::ranges(genome_gr)))
    chrom <- as.character(GenomicRanges::seqnames(genome_gr))[1]
    gene_strand <- as.character(BiocGenerics::strand(genome_gr))[1]
    toExtract_gr <- GenomicRanges::GRanges(seqnames = chrom,
        ranges = IRanges::IRanges(start - 1, width = end - start +
            3), strand = gene_strand)
    gene_gr <- IRanges::subsetByOverlaps(genome_gr, toExtract_gr)
    GenomeInfoDb::seqlevelsStyle(gene_gr) <- "UCSC"
    GenomeInfoDb::seqlevels(gene_gr) <- chrom
    gene_name_idx <- which(names(GenomicRanges::elementMetadata(gene_gr)) ==
        "gene_name")
    gene_id_idx <- which(names(GenomicRanges::elementMetadata(gene_gr)) ==
        "gene_id")
    names(GenomicRanges::elementMetadata(gene_gr))[gene_id_idx] <- "ensemble_id"
    names(GenomicRanges::elementMetadata(gene_gr))[gene_name_idx] <- "gene_id"
    gene_txdb <- GenomicFeatures::makeTxDbFromGRanges(gene_gr)
    if (label.transcripts) {
        transcript.fontsize = annotation.fontsize
    }
    else {
        transcript.fontsize = 0
    }
    gtrack <- Gviz::GeneRegionTrack(gene_txdb, start = start,
        end = end, chromosome = chrom, name = geneSymbol, just.group = "above",
        transcriptAnnotation = "symbol", showId = FALSE, fontsize.group = transcript.fontsize)
    if (!is.null(peaks.annot)) {
        start.sites <- as.numeric(sub(".*:.*:(.*)-.*:.*", "\\1",
            peaks.annot))
        end.sites <- as.numeric(sub(".*:.*:.*-(.*):.*", "\\1",
            peaks.annot))
        peak.widths <- end.sites - start.sites
        peak.names <- names(peaks.annot)
        if (is.null(peak.names))
            peak.names <- peaks.annot
        atrack <- Gviz::AnnotationTrack(start = start.sites,
            width = peak.widths, chromosome = chrom, strand = "*",
            name = "Peak", group = peak.names, genome = genome,
            just.group = "above", showId = TRUE, fontsize.group = annotation.fontsize,
            rotation.title = 90)
        gtrack <- c(gtrack, atrack)
    }
    dtrack <- list()
    wig_tracks <- list()
    if (!is.null(wig_data)) {
        if (typeof(wig_data) != "S4") {
            nc <- ncol(wig_data)
            wig_data <- GenomicRanges::makeGRangesFromDataFrame(wig_data,
                keep.extra.columns = TRUE)
        }
        GenomeInfoDb::seqlevelsStyle(wig_data) <- "UCSC"
        if (!wig_same_strand) {
            toExtract_gr <- GenomicRanges::invertStrand(toExtract_gr)
        }
        dtrack_gr <- IRanges::subsetByOverlaps(wig_data, toExtract_gr)
        GenomeInfoDb::seqlevels(dtrack_gr) <- chrom
        sample_col_idx <- 1:ncol(S4Vectors::mcols(wig_data))
        for (i in sample_col_idx) {
            tmp_gr <- dtrack_gr
            S4Vectors::mcols(tmp_gr) <- S4Vectors::mcols(tmp_gr)[i]
            dtrack_name <- names(S4Vectors::mcols(tmp_gr))
            wig_tracks[[length(wig_tracks) + 1]] <- Gviz::DataTrack(tmp_gr,
                name = dtrack_name, type = "histogram", genome = genome)
        }
    }
    if (length(bamfiles) > 0) {
        if (length(bamfile.tracknames) > 0) {
            if (length(bamfile.tracknames) == length(bamfiles)) {
                names(bamfile.tracknames) <- bamfiles
            }
            else {
                warning("BAM track names does not match number of bam files passed. \n                Replacing with filenames.")
                bamfile.tracknames <- bamfiles
                names(bamfile.tracknames) <- bamfiles
            }
        }
        else {
            bamfile.tracknames <- bamfiles
            names(bamfile.tracknames) <- bamfiles
        }
        toExtract_gr <- GenomicRanges::GRanges(seqnames = chrom,
            ranges = IRanges::IRanges(start - 50, width = end -
                start + 50), strand = gene_strand)
        for (i in bamfiles) {
            bamHeader <- Rsamtools::scanBamHeader(i)
            if (length(grep(pattern = chrom, x = names(bamHeader[[i]]$targets))) ==
                0) {
                GenomeInfoDb::seqlevelsStyle(toExtract_gr) <- "NCBI"
            }
            else {
                GenomeInfoDb::seqlevelsStyle(toExtract_gr) <- "UCSC"
            }
            param <- Rsamtools::ScanBamParam(which = toExtract_gr)
            bf <- Rsamtools::BamFile(i)
            open(bf)
            chunk0 <- GenomicAlignments::readGAlignments(bf,
                param = param)
            GenomeInfoDb::seqlevelsStyle(chunk0) <- "UCSC"
            close(bf)
            idx <- which(as.character(BiocGenerics::strand(chunk0)) ==
                gene_strand)
            if (length(idx) == 0) {
                next
            }
            tmp <- GenomicRanges::coverage(chunk0[idx])
            gr <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start:end,
                width = 1), strand = gene_strand)
            S4Vectors::mcols(gr) <- as.numeric(tmp[[chrom]])[start:end]
            dtrack[[length(dtrack) + 1]] <- Gviz::DataTrack(gr,
                name = bamfile.tracknames[i], type = "histogram",
                genome = genome)
        }
    }
    if (length(wig_tracks) > 0) {
        dtrack <- c(wig_tracks, dtrack)
    }
    toPlot <- c(gtrack, dtrack)
    if (pdf_output) {
        if (output_file_name == "") {
            warning("No file name provided")
            pdf_output = FALSE
        }
        else {
            pdf(file = output_file_name, width = 24, height = 18)
        }
    }
    extra.space = round((end - start) * 0.02)
    Gviz::plotTracks(toPlot, from = start, to = end, extend.left = extra.space,
        extend.right = extra.space, chromosome = chrom, transcriptAnnotation = "transcript",
        showId = TRUE, fontsize = axis.fontsize)
    if (zoom_3UTR) {
        idx <- which(genome_gr$type == "three_prime_utr")
        start <- min(IRanges::start(IRanges::ranges(genome_gr[idx])))
        end <- max(IRanges::end(IRanges::ranges(genome_gr[idx])))
        extra.space = round((end - start) * 1.2)
        Gviz::plotTracks(toPlot, from = end-extra.space, to = end, extend.left = extra.space,
            extend.right = extra.space, chromosome = chrom, transcriptAnnotation = "transcript",
            showId = TRUE, fontsize = axis.fontsize)
    }
    if (pdf_output) {
        dev.off()
    }
}
environment(XY_PlotCoverage) <- asNamespace("Sierra")


cal_US_in_seurat <- function(seurat=seurat,nmat=nmat){
  require(trqwe)
  require(future)
  require(future.apply)
  require(Seurat)
  both_names <- intersect(rownames(seurat@meta.data),colnames(nmat))
  nmat_tmp <- nmat[,both_names]
  US_data <- as.data.frame(nmat_tmp)
  US_data <- log(US_data+1,2)
  US_data_mean <- future_apply(US_data,2,mean)
  US_data_mean <- as.data.frame(US_data_mean)
  seurat[["US_Events"]] <- US_data_mean[rownames(seurat@meta.data),]
  return(seurat)
}



XY_rawParse <- function (data, top_genes = 50, stats = "mean")
{
    cell_group <- unique(data$cell_type)
    pb <- progress::progress_bar$new(total = length(cell_group))
    pb$tick(0)
    res <- future_lapply(cell_group,function(i){
      sub_data <- data[data$cell_type == i, ]
      counts <- t(subset(sub_data, select = -cell_type))
      counts <- apply(counts, 2, function(x) {
          storage.mode(x) <- "numeric"
          x
      })
      if (stats == "mean") {
          temp <- data.frame(rowMeans(counts), i, stringsAsFactors = FALSE)
      }
      else if (stats == "median") {
          temp <- data.frame(apply(counts, 1, FUN = median),
              i, stringsAsFactors = FALSE)
      }
      else {
          print("error stats option")
      }
      temp <- temp[order(temp[, 1], decreasing = TRUE), ]
      temp <- temp[1:ceiling(nrow(temp) * top_genes/100), ]
      temp <- temp %>% tibble::rownames_to_column()
      return(temp)
      })
    res <- do.call(rbind,res)
    pb$tick()
    colnames(res) <- c("gene", "exprs", "cell_type")
    return(res)
}
environment(XY_rawParse) <- asNamespace("iTALK")




XY_DEG <- function (data, method, min_gene_expressed = 0, min_valid_cells = 0,
    contrast = NULL, q_cut = 0.05, add = TRUE, top = 50, stats = "mean",
    ...)
{
    if (method %in% c("SCDE", "monocle", "DESingle", "MAST") &&
        dim(data)[1] >= 400) {
        print("Warning: It may take a long time. You can go and brew a cup of coffee...")
    }
    if (length(unique(data$cell_type)) != 1) {
        stop("Error: please compare data with sinlge cell type")
    }
    sub_data <- subset(data, select = -cell_type)
    combination <- combn(unique(sub_data$compare_group), 2)
    res = NULL
    if (method == "Wilcox") {
        for (i in ncol(combination)) {
            if (is.null(contrast)) {
                contrast = c(combination[, i])
            }
            sub_data <- sub_data[sub_data$compare_group %in%
                combination[, i], ]
            res <- rbind(res, WilcoxTest(sub_data, min_gene_expressed,
                min_valid_cells, contrast = contrast, ...))
        }
    }
    else if (method == "DESeq2") {
        for (i in ncol(combination)) {
            if (is.null(contrast)) {
                contrast = c(combination[, i])
            }
            sub_data <- sub_data[sub_data$compare_group %in%
                combination[, i], ]
            res <- rbind(res, DESeq2Test(sub_data, min_gene_expressed,
                min_valid_cells, contrast = contrast, parallel=TRUE,...))
        }
    }
    else if (method == "SCDE") {
        for (i in ncol(combination)) {
            if (is.null(contrast)) {
                contrast = c(combination[, i])
            }
            sub_data <- sub_data[sub_data$compare_group %in%
                combination[, i], ]
            res <- rbind(res, SCDETest(sub_data, min_gene_expressed,
                min_valid_cells, contrast = contrast, ...))
        }
    }
    else if (method == "monocle") {
        for (i in ncol(combination)) {
            if (is.null(contrast)) {
                contrast = c(combination[, i])
            }
            sub_data <- sub_data[sub_data$compare_group %in%
                combination[, i], ]
            res <- rbind(res, MonocleTest(sub_data, min_gene_expressed,
                min_valid_cells, contrast = contrast, ...))
        }
    }
    else if (method == "edgeR") {
        for (i in ncol(combination)) {
            if (is.null(contrast)) {
                contrast = c(combination[, i])
            }
            sub_data <- sub_data[sub_data$compare_group %in%
                combination[, i], ]
            res <- rbind(res, edgeRTest(sub_data, min_gene_expressed,
                min_valid_cells, contrast = contrast, ...))
        }
    }
    else if (method == "DESingle") {
        for (i in ncol(combination)) {
            if (is.null(contrast)) {
                contrast = c(combination[, i])
            }
            sub_data <- sub_data[sub_data$compare_group %in%
                combination[, i], ]
            res <- rbind(res, DESingleTest(sub_data, min_gene_expressed,
                min_valid_cells, contrast = contrast, ...))
        }
    }
    else if (method == "MAST") {
        for (i in ncol(combination)) {
            if (is.null(contrast)) {
                contrast = c(combination[, i])
            }
            sub_data <- sub_data[sub_data$compare_group %in%
                combination[, i], ]
            res <- rbind(res, MASTTest(sub_data, min_gene_expressed,
                min_valid_cells, contrast = contrast, ...))
        }
    }
    else {
        stop("Error: method currently not available")
    }
    cell_type <- unique(data$cell_type)
    res <- data.frame(res, cell_type, stringsAsFactors = FALSE)
    res <- res %>% dplyr::filter(q.value < q_cut)
    if (add) {
        parsedData <- XY_rawParse(data %>% select(-compare_group),
            top = top, stats = stats) %>% select(c(gene, cell_type)) %>%
            dplyr::mutate(logFC = 1e-04, p.value = NA, q.value = NA)
        parsedData <- parsedData %>% anti_join(res, by = c(gene = "gene"))
        res <- rbind(res, parsedData)
    }
    return(res)
}
environment(XY_DEG) <- asNamespace("iTALK")

XY_TSS_EACH_GENE <- function (object, tss.positions, assay = NULL, cells = NULL,upstream=1000,downstream=1000,
    verbose = TRUE)
{
    regions <- Extend(x = tss.positions, upstream = upstream, downstream = downstream,from.midpoint = TRUE)
    on_plus <- strand(x = regions) == "+" | strand(x = regions) =="*"
    plus.strand <- regions[on_plus, ]
    minus.strand <- regions[!on_plus, ]
    cut.matrix.plus <- MultiRegionCutMatrix(regions = plus.strand,
        object = object, assay = assay, cells = cells, verbose = FALSE)
    cut.matrix.minus <- MultiRegionCutMatrix(regions = minus.strand,
        object = object, assay = assay, cells = cells, verbose = FALSE)
    if (is.null(cut.matrix.plus)){
        full.matrix <- cut.matrix.minus[, rev(x = colnames(x = cut.matrix.minus))]
    }
    if (is.null(cut.matrix.minus)){
        full.matrix <- cut.matrix.plus
    }
    if (!is.null(cut.matrix.minus) & !is.null(cut.matrix.plus)){
        full.matrix <- cut.matrix.plus + cut.matrix.minus[, rev(x = colnames(x = cut.matrix.minus))]
    }
    colnames(full.matrix) <- -upstream:downstream
    flanking.mean <- rowMeans(x = full.matrix[, c(1:100, 1901:2001)])
    flanking.mean[flanking.mean == 0] <- mean(flanking.mean)
    norm.matrix <- full.matrix/flanking.mean
    TSS.enrichment <- rowMeans(x = norm.matrix[, 501:1500])
    return(TSS.enrichment)
}
environment(XY_TSS_EACH_GENE) <- asNamespace("Signac")


XY_subset <- function(data,colnames,sel_chr){
  sel_chr <- intersect(unique(as.character(sel_chr)),unique(as.character(data[,colnames])))
  tmp <- future_lapply(sel_chr,function(x){
    new_tmp <- data[which(data[,colnames]==x),]
    return(new_tmp)
    })
  all_sub <- do.call(rbind,tmp)
  head(all_sub)
  return(all_sub)
}



XY_gseaplot2 <- function (x, geneSetID, title = "", color = "green", base_size = 11,
    rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE,label_c=c("Description", "pvalue", "p.adjust","NES"),
    ES_geom = "line")
{
    ES_geom <- match.arg(ES_geom, c("line", "dot"))
    geneList <- position <- NULL
    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID)
    }
    else {
        gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    }
    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) +
        theme(panel.grid.major = element_line(colour = "grey92"),
            panel.grid.minor = element_line(colour = "grey92"),
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
        scale_x_continuous(expand = c(0, 0))
    if (ES_geom == "line") {
        es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description),
            size = 1)
    }
    else {
        es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description),
            size = 1, data = subset(gsdata, position == 1))
    }
    p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8),
        legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
    p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2,
            unit = "cm"))
    i <- 0
    for (term in unique(gsdata$Description)) {
        idx <- which(gsdata$ymin != 0 & gsdata$Description ==
            term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin,
        ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) +
        theme_classic(base_size) + theme(legend.position = "none",
        plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.line.x = element_blank()) +
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,
        0))
    if (length(geneSetID) == 1) {
        v <- seq(1, sum(gsdata$position), length.out = 9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0)
            inv <- inv + 1
        col = c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[unique(inv)])
        d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin,
            xmax = xmax, col = col[unique(inv)])
        p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax,
            ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d,
            alpha = 0.9, inherit.aes = FALSE)
    }
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x,
        y = ~y, yend = 0), color = "grey")
    p.pos <- p.pos + ylab("Ranked list metric") + xlab("Rank in Ordered Dataset") +
        theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2,
            l = 0.2, unit = "cm"))
    if (!is.null(title) && !is.na(title) && title != "")
        p.res <- p.res + ggtitle(title)
    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values = color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        }
        else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }
    if (pvalue_table) {
        pd <- x[geneSetID, label_c]
        pd <- pd[order(pd[, 1], decreasing = FALSE), ]
        rownames(pd) <- pd$Description
        pd <- pd[, -1]
        pd <- round(pd, 4)
        tp <- tableGrob2(pd, p.res)
        p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp,
            xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x,
                0.95), ymin = quantile(p.res$data$runningScore,
                0.75), ymax = quantile(p.res$data$runningScore,
                0.9))
    }
    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(), axis.text.x = element_text())
    if (length(subplots) == 1)
        return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2,
            r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
    if (length(rel_heights) > length(subplots))
        rel_heights <- rel_heights[subplots]
    plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}
environment(XY_gseaplot2) <- asNamespace("enrichplot")


pseudo_bulk_seurat_mean <- function(seurat_obj=seurat_obj,num_split=num_split,seed.use=seed.use,prefix=prefix,slot=slot){
  set.seed(seed.use)
  require(Seurat)
  genes.use <- rownames(seurat_obj)
  cell.sets1 <- split(colnames(seurat_obj), sort(1:length(colnames(seurat_obj))%%num_split))
  profile.set1 = matrix(, nrow = length(genes.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- GetAssayData(seurat_obj, slot = slot, assay = "RNA")[genes.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) mean(x)))
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix
    }
  }
  rownames(profile.set1) <- genes.use
  colnames(profile.set1) <- paste(prefix, 1:length(cell.sets1),sep="_")
  return(profile.set1)
}


pseudo_bulk_seurat_sum <- function(seurat_obj=seurat_obj,num_split=num_split,seed.use=seed.use,prefix=prefix,slot=slot){
  set.seed(seed.use)
  require(Seurat)
  genes.use <- rownames(seurat_obj)
  cell.sets1 <- split(colnames(seurat_obj), sort(1:length(colnames(seurat_obj))%%num_split))
  profile.set1 = matrix(, nrow = length(genes.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- GetAssayData(seurat_obj, slot = slot, assay = "RNA")[genes.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix
    }
  }
  rownames(profile.set1) <- genes.use
  colnames(profile.set1) <- paste(prefix, 1:length(cell.sets1),sep="_")
  return(profile.set1)
}


XY_runConsensusRegions <- function (testRanges, method = "majority", overlap = "any")
{
    if (length(testRanges) > 1) {
        reduced <- reduce(unlist(testRanges))
        consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
        mcols(reduced) <- do.call(cbind, lapply(testRanges, function(x) (reduced %over%
            x) + 0))
        if (method == "majority") {
            reducedConsensus <- reduced[rowSums(as.data.frame(mcols(reduced))) >
                length(testRanges)/2, ]
        }
        if (method == "none") {
            reducedConsensus <- reduced
        }
        if (is.numeric(method)) {
            reducedConsensus <- reduced[rowSums(as.data.frame(mcols(reduced))) >
                method, ]
        }
        consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
        mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)),
            consensusIDs)
        return(reducedConsensus)
    }
}
environment(XY_runConsensusRegions) <- asNamespace('soGGi')


XY_getTagMatrix <- function (chr_pos_open, weightCol = NULL, windows, flip_minor_strand = TRUE)
{
    peak.gr <- GenomicRanges::makeGRangesFromDataFrame(df=chr_pos_open,keep.extra.columns = TRUE)
    if (!is(windows, "GRanges")) {
        stop("windows should be a GRanges object...")
    }
    if (length(unique(width(windows))) != 1) {
        stop("width of windows should be equal...")
    }
    if (is.null(weightCol)) {
        peak.cov <- coverage(peak.gr)
    }
    else {
        weight <- mcols(peak.gr)[[weightCol]]
        peak.cov <- coverage(peak.gr, weight = weight)
    }
    cov.len <- elementNROWS(peak.cov)
    cov.width <- GRanges(seqnames = names(cov.len), IRanges(start = rep(1,
        length(cov.len)), end = cov.len))
    windows <- subsetByOverlaps(windows, cov.width, type = "within",
        ignore.strand = TRUE)
    chr.idx <- intersect(names(peak.cov), unique(as.character(seqnames(windows))))
    peakView <- Views(peak.cov[chr.idx], as(windows, "IntegerRangesList")[chr.idx])
    tagMatrixList <- lapply(peakView, function(x) t(viewApply(x,
        as.vector)))
    tagMatrix <- do.call("rbind", tagMatrixList)
    idx.list <- split(1:length(windows), as.factor(seqnames(windows)))
    idx <- do.call("c", idx.list)
    rownames(tagMatrix) <- idx
    tagMatrix <- tagMatrix[order(idx), ]
    if (flip_minor_strand) {
        minus.idx <- which(as.character(strand(windows)) == "-")
        tagMatrix[minus.idx, ] <- tagMatrix[minus.idx, ncol(tagMatrix):1]
    }
    tagMatrix <- tagMatrix[rowSums(tagMatrix) != 0, ]
    return(tagMatrix)
}
environment(XY_getTagMatrix) <- asNamespace('ChIPseeker')



XY_DimPlot <- function (object, dims = c(1, 2), cells = NULL, cols = NULL,
    pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL,
    shape.by = NULL, order = NULL, label = FALSE, label.size = 4,
    repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26",
    sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE)
{
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    reduction <- reduction %||% DefaultDimReduc(object = object)
    cells <- cells %||% colnames(x = object)
    data <- Embeddings(object = object[[reduction]])[cells, dims]
    data <- as.data.frame(x = data)
    dims <- paste0(Key(object = object[[reduction]]), dims)
    object[["ident"]] <- Idents(object = object)
    orig.groups <- group.by
    group.by <- group.by %||% "ident"
    data[, group.by] <- object[[group.by]][cells, , drop = TRUE]
    data <- data[order(data[,group.by]),]
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
    if (!is.null(x = split.by)) {
        ncol <- 1
    }
    if (combine) {
        plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
    }
    return(plots)
}
environment(XY_DimPlot) <- asNamespace('Seurat')

XY_seuratToURD2 <- function(seurat.object,reduction_use = reduction_use) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    # Copy over data
    ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
    if(!any(dim(seurat.object@assays$RNA@counts) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$RNA@counts[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
    # Copy over metadata
    ## TO DO - grab kmeans clustering info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) {
      get.data <- as.data.frame(seurat.object@assays$RNA@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) {
      get.data <- as.data.frame(seurat.object@meta.data)
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    # Copy over var.genes
    if(length(seurat.object@assays$RNA@var.features > 0)) ds@var.genes <- seurat.object@assays$RNA@var.features
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    }
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    } else if (.hasSlot(seurat.object, "reductions")) {
      if((reduction_use %in% names(seurat.object@reductions)) && !any(dim(Loadings(seurat.object, reduction = reduction_use)) == 0)) {
        ds@pca.load <- as.data.frame(Loadings(seurat.object, reduction = reduction_use))
        if(reduction_use=="pca") {
          ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
          ds@pca.sdev <- seurat.object@reductions$pca@stdev
          ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
        } else {
          ds@pca.scores <- as.data.frame(seurat.object@reductions$harmony@cell.embeddings)
          ds@pca.sdev <- seurat.object@reductions$harmony@stdev
          ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
        }
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
}
environment(XY_seuratToURD2) <- asNamespace('URD')

XY_RunURD_DM <- function(seurat,assay = "RNA",key = "URDDM",sigma=30,visua_group="new_anno",reduction_use="pca") {
  require(URD)
  require(Seurat)
  all_merge_URD <- XY_seuratToURD2(seurat,reduction_use=reduction_use)
  all_merge_URD <- calcDM(all_merge_URD, sigma=sigma)
  plotDimArray(all_merge_URD, reduction.use = "dm", 
    dims.to.plot = 1:8, outer.title = "URD DM", label=visua_group, plot.title="", legend=T)
  DC1 <- as.data.frame(all_merge_URD@dm$DC1)
  DC2 <- as.data.frame(all_merge_URD@dm$DC2)
  DC3 <- as.data.frame(all_merge_URD@dm$DC3)
  DC4 <- as.data.frame(all_merge_URD@dm$DC4)
  DC5 <- as.data.frame(all_merge_URD@dm$DC5)
  DC6 <- as.data.frame(all_merge_URD@dm$DC6)
  DC7 <- as.data.frame(all_merge_URD@dm$DC7)
  DC8 <- as.data.frame(all_merge_URD@dm$DC8)
  DC9 <- as.data.frame(all_merge_URD@dm$DC9)
  DC10 <- as.data.frame(all_merge_URD@dm$DC10)
  RUD_POS <- do.call(cbind,list(DC1,DC2,DC3,DC4,DC5,DC6,DC7,DC8,DC9,DC10))
  colnames(RUD_POS) <- c("DC1","DC2","DC3","DC4","DC5","DC6","DC7","DC8","DC9","DC10")
  reduction_data <- CreateDimReducObject(embeddings = as.matrix(RUD_POS),assay = assay,key = key)
  seurat[["urd"]] <- reduction_data
  return(seurat)
}


XY_new_gsea <- function(hsg_egmt_t_df_sel,hsg_genelist){
  require(DOSE)
  require(enrichplot)
  hsg_egmt_t_df_sel$ID <- 1:nrow(hsg_egmt_t_df_sel)
  rownames(hsg_egmt_t_df_sel) <- hsg_egmt_t_df_sel$ID
  geneSets <- as(hsg_egmt_t_df_sel[, "ID"], "list")
  names(geneSets) <- hsg_egmt_t_df_sel[, "ID"]
  rownames(hsg_egmt_t_df_sel) <- hsg_egmt_t_df_sel$ID
  gsea_tmp <- new("gseaResult", result = hsg_egmt_t_df_sel, geneSets = geneSets, geneList = hsg_genelist,
          params = list(pvalueCutoff = 1, nPerm = 1000,
          pAdjustMethod = "BH", exponent = 1, minGSSize = 5,
          maxGSSize = 500), readable = FALSE)
  gsea_tmp@organism <- "UNKNOWN"
  gsea_tmp@setType <- "UNKNOWN"
  gsea_tmp@keytype <- "UNKNOWN"
  return(gsea_tmp)
}



XY_ridgeplot.gseaResult <- function(x, showCategory=30, fill="p.adjust", core_enrichment = TRUE) {
    if (!is(x, "gseaResult"))
        stop("currently only support gseaResult")

    ## fill <- match.arg(fill, c("pvalue", "p.adjust", "qvalue"))
    if (fill == "qvalue") {
        fill <- "qvalues"
    }
    if (!fill %in% colnames(x@result)) {
        stop("'fill' variable not available ...")
    }

    ## geom_density_ridges <- get_fun_from_pkg('ggridges', 'geom_density_ridges')

    n <- showCategory
    if (core_enrichment) {
        gs2id <- geneInCategory(x)[seq_len(n)]
    } else {
        gs2id <- x@geneSets[x$ID[seq_len(n)]]
    }

    gs2val <- lapply(gs2id, function(id) {
        res <- x@geneList[id]
        res <- res[!is.na(res)]
    })

    nn <- names(gs2val)
    i <- match(nn, x$ID)
    nn <- x$Description[i]

    j <- order(x$NES[i], decreasing=FALSE)

    len <- sapply(gs2val, length)
    gs2val.df <- data.frame(category = rep(nn, times=len),
                            color = rep(x[i, fill], times=len),
                            value = unlist(gs2val))

    colnames(gs2val.df)[2] <- fill
    gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])
    return(gs2val.df)
  }
environment(XY_ridgeplot.gseaResult) <- asNamespace('enrichplot')


XY_make_cicero_cds <- function (cds, reduced_coordinates, k = 50, summary_stats = NULL,
    size_factor_normalize = TRUE, silent = FALSE)
{
    require(Biobase)
    assertthat::assert_that(is.data.frame(reduced_coordinates) |
        is.matrix(reduced_coordinates))
    assertthat::assert_that(assertthat::are_equal(nrow(reduced_coordinates),
        nrow(pData(cds))))
    assertthat::assert_that(setequal(row.names(reduced_coordinates),
        colnames(cds)))
    assertthat::assert_that(assertthat::is.count(k) & k > 1)
    assertthat::assert_that(is.character(summary_stats) | is.null(summary_stats))
    if (!is.null(summary_stats)) {
        assertthat::assert_that(all(summary_stats %in% names(pData(cds))),
            msg = paste("One of your summary_stats is missing",
                "from your pData table. Either add a", "column with the name in",
                "summary_stats, or remove the name", "from the summary_stats parameter.",
                collapse = " "))
        assertthat::assert_that(sum(vapply(summary_stats, function(x) {
            !(is(pData(cds)[, x], "numeric") | is(pData(cds)[,
                x], "integer"))
        }, 1)) == 0, msg = paste("All columns in summary_stats must be",
            "of class numeric or integer.", collapse = " "))
    }
    assertthat::assert_that(is.logical(size_factor_normalize))
    assertthat::assert_that(is.logical(silent))
    reduced_coordinates <- as.data.frame(reduced_coordinates)
    reduced_coordinates <- reduced_coordinates[colnames(cds),
        ]
    nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
    row.names(nn_map) <- row.names(reduced_coordinates)
    nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
    good_choices <- seq_len(nrow(nn_map))
    choice <- sample(seq_len(length(good_choices)), size = 1,
        replace = FALSE)
    chosen <- good_choices[choice]
    good_choices <- good_choices[good_choices != good_choices[choice]]
    it <- 0
    k2 <- k * 2
    get_shared <- function(other, this_choice) {
        k2 - length(union(cell_sample[other, ], this_choice))
    }
    while (length(good_choices) > 0 & it < 5000) {
        it <- it + 1
        choice <- sample(seq_len(length(good_choices)), size = 1,
            replace = FALSE)
        new_chosen <- c(chosen, good_choices[choice])
        good_choices <- good_choices[good_choices != good_choices[choice]]
        cell_sample <- nn_map[new_chosen, ]
        others <- seq_len(nrow(cell_sample) - 1)
        this_choice <- cell_sample[nrow(cell_sample), ]
        shared <- sapply(others, get_shared, this_choice = this_choice)
        if (max(shared) < 0.9 * k) {
            chosen <- new_chosen
        }
    }
    cell_sample <- nn_map[chosen, ]
    if (!silent) {
        combs <- combn(nrow(cell_sample), 2)
        shared <- apply(combs, 2, function(x) {
            k2 - length(unique(as.vector(cell_sample[x, ])))
        })
        message(paste0("Overlap QC metrics:\nCells per bin: ",
            k, "\nMaximum shared cells bin-bin: ", max(shared),
            "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
            median(shared)))
        if (mean(shared)/k > 0.1)
            warning("On average, more than 10% of cells are shared between paired bins.")
    }
    exprs_old <- exprs(cds)
    mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
        cell_sample[x, , drop = FALSE])
    mask <- Matrix::Matrix(mask)
    new_exprs <- exprs_old %*% mask
    new_exprs <- Matrix::t(new_exprs)
    new_exprs <- as.matrix(new_exprs)
    pdata <- pData(cds)
    new_pcols <- "agg_cell"
    if (!is.null(summary_stats)) {
        new_pcols <- c(new_pcols, paste0("mean_", summary_stats))
    }
    new_pdata <- plyr::adply(cell_sample, 1, function(x) {
        sub <- pdata[x, ]
        df_l <- list()
        df_l["temp"] <- 1
        for (att in summary_stats) {
            df_l[paste0("mean_", att)] <- mean(sub[, att])
        }
        data.frame(df_l)
    })
    new_pdata$agg_cell <- paste("agg", chosen, sep = "")
    new_pdata <- new_pdata[, new_pcols, drop = FALSE]
    row.names(new_pdata) <- new_pdata$agg_cell
    row.names(new_exprs) <- new_pdata$agg_cell
    new_exprs <- as.matrix(t(new_exprs))
    fdf <- fData(cds)
    new_pdata$temp <- NULL
    fd <- new("AnnotatedDataFrame", data = as.data.frame(fdf))
    pd <- new("AnnotatedDataFrame", data = new_pdata)
    cicero_cds <- suppressWarnings(newCellDataSet(new_exprs,
        phenoData = pd, featureData = fd, expressionFamily = negbinomial.size(),
        lowerDetectionLimit = 0))
    cicero_cds <- monocle::detectGenes(cicero_cds, min_expr = 0.1)
    cicero_cds <- BiocGenerics::estimateSizeFactors(cicero_cds)
    if (any(!c("chr", "bp1", "bp2") %in% names(cicero_cds@featureData@data))) {
        cicero_cds@featureData@data$chr <- NULL
        cicero_cds@featureData@data$bp1 <- NULL
        cicero_cds@featureData@data$bp2 <- NULL
        cicero_cds@featureData@data <- cbind(cicero_cds@featureData@data, df_for_coords(row.names(cicero_cds@featureData@data)))
    }
    if (size_factor_normalize) {
        Biobase::exprs(cicero_cds) <- t(t(Biobase::exprs(cicero_cds))/Biobase::pData(cicero_cds)$Size_Factor)
    }
    cicero_cds
}
environment(XY_make_cicero_cds) <- asNamespace('monocle3')

XY_plotPCA = function(data, meta_info=meta_info,intgroup="condition", ntop=500, returnData=FALSE,label=label)
{
  # calculate the variance for each gene
  rv <- rowVars(data)
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(data[select,])
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  intgroup.df <- as.data.frame(meta_info[, intgroup, drop=FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    meta_info[[intgroup]]
  }
  d <- data.frame(PC1=pca$rotation[,1], PC2=pca$rotation[,2], group=group, intgroup.df, name=colnames(data))
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +theme(legend.position="none")+
      geom_text_repel(
        data = d,
        aes(label = d[,label]),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )
}
environment(XY_plotPCA) <- asNamespace('DESeq2')



require(Seurat)
require(patchwork)
require(ggplot2)
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot <- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



XY_DUTest <- function (apa.seurat.object, population.1, population.2 = NULL,
    exp.thresh = 0.1, fc.thresh = 0.25, adj.pval.thresh = 0.05,
    num.splits = 6, seed.use = 1, feature.type = c("UTR3", "UTR5",
        "exon", "intron"), include.annotations = FALSE, filter.pA.stretch = FALSE,
    verbose = TRUE, do.MAPlot = FALSE, return.dexseq.res = FALSE,
    ncores = 1)
{
    if (!"DEXSeq" %in% rownames(x = installed.packages())) {
        stop("Please install DEXSeq before using this function\n         (http://bioconductor.org/packages/release/bioc/html/DEXSeq.html)")
    }
    high.expressed.peaks <- GetExpressedPeaks(apa.seurat.object,
        population.1, population.2, threshold = exp.thresh)
    length(high.expressed.peaks)
    annot.subset <- Tool(apa.seurat.object, "Sierra")[high.expressed.peaks,
        ]
    peaks.to.use <- apply(annot.subset, 1, function(x) {
        ifelse(sum(x[feature.type] == "YES") >= 1, TRUE, FALSE)
    })
    peaks.to.use <- names(peaks.to.use[which(peaks.to.use ==
        TRUE)])
    high.expressed.peaks <- intersect(high.expressed.peaks, peaks.to.use)
    if (verbose)
        print(paste(length(high.expressed.peaks), "expressed peaks in feature types",
            toString(feature.type)))
    if (filter.pA.stretch) {
        if (is.null(Tool(apa.seurat.object, "Sierra")$pA_stretch)) {
            stop("pA_stretch not in annotation data: please run nnotate_gr_from_gtf with\n           an input genome to provide required annotation.")
        } else {
            annot.subset <- Tool(apa.seurat.object, "Sierra")[high.expressed.peaks,
                ]
            peaks.non.arich <- rownames(subset(annot.subset,
                pA_stretch == FALSE))
            high.expressed.peaks <- intersect(high.expressed.peaks,
                peaks.non.arich)
            if (verbose)
                print(paste(length(high.expressed.peaks), "peaks after filtering out A-rich annotations"))
        }
    }
    annotations.highly.expressed <- Tool(apa.seurat.object, "Sierra")[high.expressed.peaks,
        ]
    annotations.highly.expressed <- subset(annotations.highly.expressed,
        Gene_name != "")
    high.expressed.peaks <- rownames(annotations.highly.expressed)
    gene.names <- annotations.highly.expressed[, "Gene_name"]
    gene.table <- table(gene.names)
    multi.genes <- gene.table[gene.table > 1]
    if (verbose)
        print(paste(length(multi.genes), "genes detected with multiple peak sites expressed"))
    multi.gene.names <- names(multi.genes)
    peaks.use <- high.expressed.peaks[which(gene.names %in% multi.gene.names)]
    if (verbose)
        print(paste(length(peaks.use), "individual peak sites to test"))
    set.seed(seed.use)
    if (length(population.1) == 1) {
        cells.1 <- names(Seurat::Idents(apa.seurat.object))[which(Seurat::Idents(apa.seurat.object) ==
            population.1)]
    } else {
        cells.1 <- population.1
    }
    cells.1 = sample(cells.1)
    cell.sets1 <- split(cells.1, sort(1:length(cells.1)%%num.splits))
    profile.set1 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets1))
    for (i in 1:length(cell.sets1)) {
        this.set <- cell.sets1[[i]]
        sub.matrix <- Seurat::GetAssayData(apa.seurat.object,
            slot = "counts", assay = "RNA")[peaks.use, this.set]
        if (length(this.set) > 1) {
            this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
            profile.set1[, i] <- this.profile
        } else {
            profile.set1[, i] <- sub.matrix
        }
    }
    rownames(profile.set1) <- peaks.use
    colnames(profile.set1) <- paste0("Population1_", 1:length(cell.sets1))
    if (is.null(population.2)) {
        cells.2 <- setdiff(colnames(apa.seurat.object), cells.1)
    } else {
        if (length(population.2) == 1) {
            cells.2 <- names(Seurat::Idents(apa.seurat.object))[which(Seurat::Idents(apa.seurat.object) ==
                population.2)]
        } else {
            cells.2 <- population.2
        }
    }
    cells.2 = sample(cells.2)
    cell.sets2 <- split(cells.2, sort(1:length(cells.2)%%num.splits))
    profile.set2 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets2))
    for (i in 1:length(cell.sets2)) {
        this.set <- cell.sets2[[i]]
        sub.matrix <- Seurat::GetAssayData(apa.seurat.object,
            slot = "counts", assay = "RNA")[peaks.use, this.set]
        if (length(this.set) > 1) {
            this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
            profile.set2[, i] <- this.profile
        } else {
            profile.set2[, i] <- sub.matrix
        }
    }
    rownames(profile.set2) <- peaks.use
    colnames(profile.set2) <- paste0("Population2_", 1:length(cell.sets2))
    peak.matrix <- cbind(profile.set1, profile.set2)
    sampleTable <- data.frame(row.names = c(colnames(profile.set1),
        colnames(profile.set2)), condition = c(rep("target",
        ncol(profile.set1)), rep("comparison", ncol(profile.set2))))
    dexseq.feature.table <- Tool(apa.seurat.object, "Sierra")[,
        c("Gene_name", "Gene_part", "Peak_number")]
    dexseq.feature.table$Peak <- rownames(dexseq.feature.table)
    dexseq.feature.table <- dexseq.feature.table[rownames(peak.matrix),
        ]
    rownames(dexseq.feature.table) <- paste0(dexseq.feature.table$Gene_name,
        ":", dexseq.feature.table$Peak_number)
    rownames(peak.matrix) <- rownames(dexseq.feature.table)
    peak_ID_set = dexseq.feature.table[rownames(peak.matrix),
        "Peak_number"]
    gene_names = dexseq.feature.table[rownames(peak.matrix),
        "Gene_name"]
    dxd = DEXSeq::DEXSeqDataSet(peak.matrix, sampleData = sampleTable,
        groupID = gene_names, featureID = peak_ID_set, design = ~sample +
            exon + condition:exon)
    if (verbose)
        print("Running DEXSeq test...")
    if (ncores > 1) {
        BPPARAM = BiocParallel::MulticoreParam(workers = ncores)
        dxd = DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
        dxd = DEXSeq::estimateDispersions(dxd, BPPARAM = BPPARAM)
        dxd = DEXSeq::testForDEU(dxd, BPPARAM = BPPARAM)
        dxd = DEXSeq::estimateExonFoldChanges(dxd, BPPARAM = BPPARAM)
        dxr1 = DEXSeq::DEXSeqResults(dxd)
    } else {
        dxd = DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
        dxd = DEXSeq::estimateDispersions(dxd)
        dxd = DEXSeq::testForDEU(dxd)
        dxd = DEXSeq::estimateExonFoldChanges(dxd)
        dxr1 = DEXSeq::DEXSeqResults(dxd)
    }
    if (do.MAPlot)
        DEXSeq::plotMA(dxr1, alpha = adj.pval.thresh, ylim = c(min(dxr1$log2fold_target_comparison),
            max(dxr1$log2fold_target_comparison)))
    if (return.dexseq.res)
        return(dxr1)
    dxrSig <- subset(as.data.frame(dxr1), padj < adj.pval.thresh &
        abs(log2fold_target_comparison) > fc.thresh)
    dxrSig_subset <- dxrSig[, c("groupID", "exonBaseMean", "pvalue",
        "padj", "log2fold_target_comparison")]
    peaks.to.add = dexseq.feature.table[rownames(dxrSig_subset),
        "Peak"]
    rownames(dxrSig_subset) = peaks.to.add
    population.1.pct <- get_percent_expression(apa.seurat.object,
        population.1, remainder = FALSE, geneSet = rownames(dxrSig_subset))
    if (is.null(population.2)) {
        population.2.pct <- get_percent_expression(apa.seurat.object,
            population.1, remainder = TRUE, geneSet = rownames(dxrSig_subset))
        population.2 <- "Remainder"
    } else {
        population.2.pct <- get_percent_expression(apa.seurat.object,
            population.2, remainder = FALSE, geneSet = rownames(dxrSig_subset))
    }
    dxrSig_subset$population1_pct <- population.1.pct
    dxrSig_subset$population2_pct <- population.2.pct
    feature.type <- Tool(apa.seurat.object, "Sierra")[rownames(dxrSig_subset),
        c("FeaturesCollapsed")]
    dxrSig_subset$feature_type = feature.type
    if (include.annotations) {
        junction.annot <- Tool(apa.seurat.object, "Sierra")[rownames(dxrSig_subset),
            c("pA_motif", "pA_stretch")]
        dxrSig_subset <- cbind(dxrSig_subset, junction.annot)
        dxrSig_subset <- dxrSig_subset[, c("groupID", "feature_type",
            "pA_motif", "pA_stretch", "population1_pct",
            "population2_pct", "pvalue", "padj", "log2fold_target_comparison")]
        colnames(dxrSig_subset) <- c("gene_name", "genomic_feature(s)",
            "pA_motif", "pA_stretch", "population1_pct",
            "population2_pct", "pvalue", "padj", "Log2_fold_change")
    } else {
        dxrSig_subset <- dxrSig_subset[, c("groupID", "feature_type",
            "population1_pct", "population2_pct", "pvalue", "padj",
            "log2fold_target_comparison")]
        colnames(dxrSig_subset) <- c("gene_name", "genomic_feature(s)",
            "population1_pct", "population2_pct", "pvalue", "padj",
            "Log2_fold_change")
    }
    dxrSig_subset <- dxrSig_subset[order(dxrSig_subset$padj,
        decreasing = FALSE), ]
    return(dxrSig_subset)
}
environment(XY_DUTest) <- asNamespace('Sierra')



XY_PlotCoverage_selected <- function (genome_gr, geneSymbol = "",chrom,start,end, bamfiles = NULL, peaks.annot = NULL,
 genome = NULL,annotation.fontsize = NULL) {
    GenomeInfoDb::seqlevelsStyle(genome_gr) <- "UCSC"
    idx <- which(genome_gr$gene_name == geneSymbol)
    genome_gr <- genome_gr[idx]
    start_sel <- min(IRanges::start(IRanges::ranges(genome_gr)))
    end_sel <- max(IRanges::end(IRanges::ranges(genome_gr)))
    chrom_sel <- as.character(GenomicRanges::seqnames(genome_gr))[1]
    gene_strand <- as.character(BiocGenerics::strand(genome_gr))[1]
    toExtract_gr <- GenomicRanges::GRanges(seqnames = chrom_sel,ranges = IRanges::IRanges(start_sel - 1, width = end_sel - start_sel +3), strand = gene_strand)
    gene_gr <- IRanges::subsetByOverlaps(genome_gr, toExtract_gr)
    sel_gr <- GenomicRanges::GRanges(seqnames = chrom,ranges = IRanges::IRanges(start - 1, width = end - start +3), strand = gene_strand)
    tmp_gr <- findOverlaps(sel_gr, gene_gr, type="any")
    gene_gr <- gene_gr[tmp_gr@to,]
    GenomeInfoDb::seqlevelsStyle(gene_gr) <- "UCSC"
    GenomeInfoDb::seqlevels(gene_gr) <- chrom
    gene_name_idx <- which(names(GenomicRanges::elementMetadata(gene_gr)) =="gene_name")
    gene_id_idx <- which(names(GenomicRanges::elementMetadata(gene_gr)) =="gene_id")
    names(GenomicRanges::elementMetadata(gene_gr))[gene_id_idx] <- "ensemble_id"
    names(GenomicRanges::elementMetadata(gene_gr))[gene_name_idx] <- "gene_id"
    gene_txdb <- GenomicFeatures::makeTxDbFromGRanges(gene_gr)
    transcript.fontsize = annotation.fontsize
    gtrack <- Gviz::GeneRegionTrack(gene_txdb, start = start,
        end = end, chromosome = chrom, name = geneSymbol, just.group = "above",
        transcriptAnnotation = "symbol", showId = FALSE, fontsize.group = transcript.fontsize)
    start.sites <- as.numeric(sub(".*:.*:(.*)-.*:.*", "\\1",peaks.annot))
    end.sites <- as.numeric(sub(".*:.*:.*-(.*):.*", "\\1",peaks.annot))
    peak.widths <- end.sites - start.sites
    peak.names <- names(peaks.annot)
    peak.names <- peaks.annot
    atrack <- Gviz::AnnotationTrack(start = start.sites,
        width = peak.widths, chromosome = chrom, strand = "*",
        name = "Peak", group = peak.names, genome = genome,
        just.group = "above", showId = TRUE, fontsize.group = annotation.fontsize,
        rotation.title = 90)
    gtrack <- c(gtrack, atrack)
    dtrack <- list()
    bamfile.tracknames <- bamfiles
    names(bamfile.tracknames) <- bamfiles
    toExtract_gr <- GenomicRanges::GRanges(seqnames = chrom,ranges = IRanges::IRanges(start - 50, width = end - start + 50), strand = gene_strand)
    for (i in bamfiles) {
        bamHeader <- Rsamtools::scanBamHeader(i)
        if (length(grep(pattern = chrom, x = names(bamHeader[[i]]$targets))) ==0) {
            GenomeInfoDb::seqlevelsStyle(toExtract_gr) <- "NCBI"
        } else {
            GenomeInfoDb::seqlevelsStyle(toExtract_gr) <- "UCSC"
        }
        param <- Rsamtools::ScanBamParam(which = toExtract_gr)
        bf <- Rsamtools::BamFile(i)
        open(bf)
        chunk0 <- GenomicAlignments::readGAlignments(bf,param = param)
        GenomeInfoDb::seqlevelsStyle(chunk0) <- "UCSC"
        close(bf)
        idx <- which(as.character(BiocGenerics::strand(chunk0)) == gene_strand)
        if (length(idx) == 0) {
            next
        }
        tmp <- GenomicRanges::coverage(chunk0[idx])
        gr <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start:end,width = 1), strand = gene_strand)
        S4Vectors::mcols(gr) <- as.numeric(tmp[[chrom]])[start:end]
        dtrack[[length(dtrack) + 1]] <- Gviz::DataTrack(gr,name = bamfile.tracknames[i], type = "histogram",genome = genome)
    }
    toPlot <- c(gtrack, dtrack)
    extra.space = round((end - start) * 0.02)
    Gviz::plotTracks(toPlot, from = start, to = end, extend.left = extra.space,
        extend.right = extra.space, chromosome = chrom, transcriptAnnotation = "transcript",
        showId = TRUE, fontsize = NULL)
}
environment(XY_PlotCoverage) <- asNamespace("Sierra")

Binner <- function(cds_object,cells_subset,anno_group){
  df <- data.frame(pData(cds_object[,cells_subset]))
  df <- df[,c("Pseudotime", anno_group)]
  colnames(df) <- c("Pseudotime", "State")
  df <- df[order(df$Pseudotime, decreasing = F),]
  len <- length(df$Pseudotime)
  bin <- round(len/100)
  State <- c()
  value <- c()
  for(i in 0:99){
    if(i < 99){
      start <- 1+(bin*i)
      stop <- bin+(bin*i)
      value <- df$State[c(start:stop)][length(df$State[c(start:stop)])/2]
      State <- c(State, value)
    }
    else{
      State <- c(State, value)
    }
  }
  return(as.data.frame(State))
}


enrich_pvalue <- function(N, A, B, k)
{
  require(gmp)
    m <- A + k
    n <- B + k
    i <- k:min(m,n)

    as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}

XY_LRPlot <- function (data, datatype, gene_col = NULL, transparency = 0.5,
    link.arr.lwd = 1, link.arr.lty = NULL, link.arr.col = NULL,
    link.arr.width = NULL, link.arr.type = NULL, facing = "clockwise",
    cell_col = NULL, print.cell = TRUE, track.height_1 = uh(2,
        "mm"), track.height_2 = uh(12, "mm"), annotation.height_1 = 0.01,
    annotation.height_2 = 0.01, text.vjust = "0.4cm", ...)
{
    cell_group <- unique(c(data$cell_from, data$cell_to))
    genes <- c(structure(data$ligand, names = data$cell_from),
        structure(data$receptor, names = data$cell_to))
    genes <- genes[!duplicated(paste(names(genes), genes))]
    genes <- genes[order(names(genes))]
    if (is.null(link.arr.lty)) {
        if (datatype == "mean count") {
            link.arr.lty = "solid"
        }
        else if (datatype == "DEG") {
            link.arr.lty = structure(ifelse(data$cell_from_logFC ==
                1e-04, "dashed", "solid"), names = paste(data$cell_from,
                data$receptor))
        }
        else {
            print("invalid datatype")
        }
    }
    if (is.null(link.arr.col)) {
        if (datatype == "mean count") {
            data <- data %>% mutate(link_col = "black")
        }
        else if (datatype == "DEG") {
            data <- data %>% mutate(link_col = ifelse(cell_from_logFC ==
                1e-04, ifelse(cell_to_logFC > 0, "#d73027", "#d73027"),
                ifelse(cell_to_logFC == 1e-04, ifelse(cell_from_logFC >
                  0, "#d73027", "#d73027"), ifelse(cell_from_logFC >
                  0, ifelse(cell_to_logFC > 0, "#d73027", "#d73027"),
                  ifelse(cell_to_logFC > 0, "#d73027", "#d73027")))))
            data[data$cell_to_logFC==0,]$link_col <- "#ffffff"
            data[data$cell_from_logFC==0 ,]$link_col <- "#ffffff"
        }
        else {
            print("invalid datatype")
        }
    }
    else {
        data$link_col = link.arr.col
    }
    if (is.null(link.arr.type)) {
        if (datatype == "mean count") {
            link.arr.type = "triangle"
        }
        else if (datatype == "DEG") {
            link.arr.type = structure(ifelse(data$cell_to_logFC ==
                1e-04, "ellipse", "triangle"), names = paste(data$cell_from,
                data$receptor))
        }
        else {
            print("invalid datatype")
        }
    }
    if (is.null(gene_col)) {
        comm_col <- structure(c("#99ff99", "#99ccff", "#ff9999",
            "#ffcc99"), names = c("other", "cytokine", "checkpoint",
            "growth factor"))
        gene_col <- structure(c(comm_col[data$comm_type], rep("#073c53",
            length(data$receptor))), names = c(data$ligand, data$receptor))
    }
    if (is.null(cell_col)) {
        cell_col <- structure(randomColor(count = length(unique(names(genes))),
            luminosity = "dark"), names = unique(names(genes)))
    }
    if (is.null(link.arr.lwd)) {
        data <- data %>% mutate(arr_width = 1)
        data[data$cell_to_logFC==0,]$arr_width <- 0
        data[data$cell_from_logFC==0,]$arr_width <- 0
    }
    else if (max(abs(link.arr.lwd)) - min(abs(link.arr.lwd)) ==
        0 && all(link.arr.lwd != 1e-04)) {
        data <- data %>% mutate(arr_width = ifelse(abs(link.arr.lwd <
            5), abs(link.arr.lwd), 5))
        data[data$cell_to_logFC==0,]$arr_width <- 0
        data[data$cell_from_logFC==0,]$arr_width <- 0
    }
    else {
        data <- data %>% mutate(arr_width = ifelse(link.arr.lwd ==
            1e-04, 2, 1 + 5/(max(abs(link.arr.lwd)) - min(abs(link.arr.lwd))) *
            (abs(link.arr.lwd) - min(abs(link.arr.lwd)))))
        data[data$cell_to_logFC==0,]$arr_width <- 0
        data[data$cell_from_logFC==0,]$arr_width <- 0
    }
    if (length(cell_group) != 1) {
        gap.degree <- do.call("c", lapply(table(names(genes)),
            function(i) c(rep(1, i - 1), 8)))
    }
    else {
        gap.degree <- do.call("c", lapply(table(names(genes)),
            function(i) c(rep(1, i))))
    }
    circos.par(gap.degree = gap.degree)
    if (length(gene_col) == 1) {
        grid.col = gene_col
    }
    else {
        grid.col = gene_col[genes]
        names(grid.col) <- paste(names(genes), genes)
    }
    if (is.null(link.arr.width)) {
        data <- data %>% mutate(link.arr.width = data$arr_width/10)
    }
    else if (max(abs(link.arr.width)) - min(abs(link.arr.width)) ==
        0 && all(link.arr.width != 1e-04)) {
        data <- data %>% mutate(link.arr.width = ifelse(abs(link.arr.width) <
            0.5, abs(link.arr.width), 0.5))
        data[data$cell_to_logFC==0,]$link.arr.width <- 0
        data[data$cell_from_logFC==0,]$link.arr.width <- 0
    }
    else {
        data <- data %>% mutate(link.arr.width = ifelse(link.arr.width ==
            1e-04, 0.2, (1 + 5/(max(abs(link.arr.width)) - min(abs(link.arr.width))) *
            (abs(link.arr.width) - min(abs(link.arr.width))))/10))
        data[data$cell_to_logFC==0,]$link.arr.width <- 0
        data[data$cell_from_logFC==0,]$link.arr.width <- 0
    }
    data$link_col <- factor(data$link_col,levels=c("#d73027","#ffffff"))
    data <- data[order(data$link_col,decreasing=TRUE),]
    data$link_col <- as.character(data$link_col)
    chordDiagram(as.data.frame(cbind(paste(data$cell_from, data$ligand),
        paste(data$cell_to, data$receptor))), order = paste(names(genes),
        genes), grid.col = grid.col, transparency = transparency,
        directional = 1, direction.type = "arrows", link.arr.lwd = data$arr_width,
        link.arr.lty = link.arr.lty, link.arr.type = link.arr.type,
        link.arr.width = data$link.arr.width, link.arr.col = data$link_col,
        col = "#00000000", annotationTrack = c("grid"), preAllocateTracks = list(list(track.height = track.height_1),
            list(track.height = track.height_2)), annotationTrackHeight = c(annotation.height_1,
            annotation.height_2), ...)
    circos.trackPlotRegion(track.index = 2, panel.fun = function(x,
        y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.index = genes[get.cell.meta.data("sector.numeric.index")]
        circos.text(mean(xlim), mean(ylim), sector.index, col = "black",
            cex = 0.7, facing = facing, niceFacing = TRUE)
    }, bg.border = 0)
    if (print.cell) {
        for (c in unique(names(genes))) {
            gene = as.character(genes[names(genes) == c])
            highlight.sector(sector.index = paste(c, gene), track.index = 1,
                col = ifelse(length(cell_col) == 1, cell_col,
                  cell_col[c]), text = c, text.vjust = text.vjust,
                niceFacing = TRUE, lwd = 1)
        }
    }
    circos.clear()
}
environment(XY_LRPlot) <- asNamespace("iTALK")


miQC.keep_info <- function (sce, model = NULL, posterior_cutoff = 0.75, keep_all_below_boundary = TRUE,
    enforce_left_cutoff = TRUE, verbose = TRUE)
{
    metrics <- as.data.frame(colData(sce))
    if (is.null(model)) {
        warning("call 'mixtureModel' explicitly to get stable model features")
        model <- mixtureModel(sce)
    }
    intercept1 <- parameters(model, component = 1)[1]
    intercept2 <- parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
        compromised_dist <- 1
        intact_dist <- 2
    }
    else {
        intact_dist <- 1
        compromised_dist <- 2
    }
    post <- posterior(model)
    metrics$prob_compromised <- post[, compromised_dist]
    sce$prob_compromised <- metrics$prob_compromised
    metrics$keep <- metrics$prob_compromised <= posterior_cutoff
    if (keep_all_below_boundary == TRUE) {
        predictions <- fitted(model)[, intact_dist]
        metrics$intact_prediction <- predictions
        metrics[metrics$subsets_mito_percent < metrics$intact_prediction,
            ]$keep <- TRUE
    }
    if (enforce_left_cutoff == TRUE) {
        min_discard <- min(metrics[!metrics$keep, ]$subsets_mito_percent)
        min_index <- which(metrics$subsets_mito_percent == min_discard)[1]
        lib_complexity <- metrics[min_index, ]$detected
        metrics[metrics$detected <= lib_complexity & metrics$subsets_mito_percent >=
            min_discard, ]$keep <- FALSE
    }
    if (verbose == TRUE) {
        to_remove <- length(which(metrics$keep == FALSE))
        total <- length(metrics$keep)
        cat("Removing", to_remove, "out of", total, "cells.")
    }
    return(metrics)
}
environment(miQC.keep_info) <- asNamespace('miQC')


XY_scale <- function(obj=obj,min_max_num=i) {
  chonglai_zscore_1 <- t(apply(obj, 1, function(x) (x-mean(x))/sd(x)))
  # chonglai_zscore_1 <- na.omit(chonglai_zscore_1)
  test <- as.data.frame(chonglai_zscore_1)
  test[test > min_max_num] <- min_max_num
  test[test < -min_max_num] <- -min_max_num
  return(test)
}

venntest_three <- function(N,m,a,b,c){
#N=7000 # total number of genes for selection
#m=3 # number of overlaps in all three sets
#a=200 # number of up-regulated genes in sample A
#b=250 # number of up-regulated genes in sample B
#c=300 # number of up-regualted genes in sample C
#create a dataframe to store frequency table
d<-data.frame(0:min(a,b,c),rep(0,min(a,b,c)+1))
for (i in 1:10000){
A=sample(1:N,size=a,replace=FALSE)
B=sample(1:N,size=b,replace=FALSE)
C=sample(1:N,size=c,replace=FALSE)
e<-table((C %in% A)&(C %in% B))["TRUE"]
# if there is no intersection, put 0 instead of NA
if(!complete.cases(e)){
e<-0
}
# add one to the counter for corresponding occurence
d[e+1,2]<-d[e+1,2]+1
}
colnames(d)<-c("Intersect","p-value")
# convert counts to ratios
d[,2]<-d[,2]/10000
# calculate p-value
p<-sum(d[(m+1):(min(a,b,c)+1),2])
#print(d)
return(p)
}


XY_covplot <- function (peak, weightCol = NULL, xlab = "Chromosome Size (bp)",
    ylab = "", title = "ChIP Peaks over Chromosomes", chrs = NULL,
    xlim = NULL, lower = 1) {
    isList <- FALSE
    if (is(peak, "GRanges") || length(peak) == 1) {
        tm <- ChIPseeker:::getChrCov(peak = peak, weightCol = weightCol, chrs,
            xlim, lower = lower)
    } else {
        isList <- TRUE
        ltm <- lapply(peak, ChIPseeker:::getChrCov, weightCol = weightCol,
            chrs = chrs, xlim = xlim, lower = lower)
        if (is.null(names(ltm))) {
            nn <- paste0("peak", seq_along(ltm))
            warning("input is not a named list, set the name automatically to ",
                paste(nn, collapse = " "))
            names(ltm) <- nn
        }
        tm <- ChIPseeker:::list_to_dataframe(ltm)
        chr.sorted <- ChIPseeker:::sortChrName(as.character(unique(tm$chr)))
        tm$chr <- factor(tm$chr, levels = chr.sorted)
    }
    chr <- start <- end <- value <- .id <- NULL
    tm$value <- log(tm$value)   
    p <- ggplot(tm, aes(start, value))
    if (isList) {
        p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end,
            ymax = value, fill = .id, color = .id))
    }
    else {
        p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end,
            ymax = value), fill = "black", color = "black")
    }
    if (length(unique(tm$chr)) > 1) {
        p <- p + facet_grid(chr ~ ., scales = "free")
    }
    p <- p + theme_classic()
    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
    p <- p + scale_y_continuous(expand = c(0, 0))
    p <- p + theme(strip.text.y = element_text(angle = 360))
    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) &&
        length(xlim) == 2) {
        p <- p + xlim(xlim)
    }
    return(p)
}
environment(XY_covplot) <- asNamespace('ChIPseeker')

XY_convert_DF_to_GRange <- function(peak.df) {
  if (unique(unique(peak.df[, 5]) %in% c("+", "-", "*"))==TRUE ){
    peak.gr = GRanges(seqnames = peak.df[, 1], ranges = IRanges(peak.df[,2], peak.df[, 3]))
  } else {
    peak.gr = GRanges(seqnames = peak.df[, 1], ranges = IRanges(peak.df[,2], peak.df[, 3]), strand ="*")
  }
  cn <- colnames(peak.df)
  if (length(cn) > 3) {
      for (i in 4:length(cn)) {
          mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
      }
  }
  return(peak.gr)
}

XY_circos.genomicDensity <- function (data, y_range=y_range, window.size = NULL, overlap = TRUE,
    count_by = c("percent", "number"), col = ifelse(area, "grey",
        "black"), lwd = par("lwd"), lty = par("lty"), type = "l",
    area = TRUE, area.baseline = NULL, baseline = 0, border = NA,
    ...)
{
    if (!is.null(area.baseline)) {
        baseline = area.baseline
        warning_wrap("`area.baseline` is deprecated, please use `baseline` instead.")
    }
    data = normalizeToDataFrame(data)
    if (!is.dataFrameList(data)) {
        data = list(data)
    }
    if (length(col) == 1) {
        col = rep(col, length(data))
    }
    if (length(lwd) == 1) {
        lwd = rep(lwd, length(data))
    }
    if (length(lty) == 1) {
        lty = rep(lty, length(data))
    }
    if (length(type) == 1) {
        type = rep(type, length(data))
    }
    if (length(area) == 1) {
        area = rep(area, length(data))
    }
    if (length(baseline) == 1) {
        baseline = rep(baseline, length(data))
    }
    if (length(border) == 1) {
        border = rep(border, length(data))
    }
    s = sapply(get.all.sector.index(), function(si) get.cell.meta.data("xrange",
        sector.index = si))
    if (is.null(window.size)) {
        window.size = 10^nchar(sum(s))/1000
    }
    df = vector("list", length = length(data))
    for (i in seq_along(data)) {
        df[[i]] = genomicDensity(data[[i]], window.size = window.size,
            overlap = overlap, count_by = count_by)
    }
    if (length(df) == 1) {
        circos.genomicTrackPlotRegion(df[[1]], ylim = c(min(max(y_range)), max(y_range)),
            panel.fun = function(region, value, ...) {
                circos.genomicLines(region, value, col = col,
                  lwd = lwd, lty = lty, type = type, border = border,
                  area = area, baseline = baseline)
            }, ...)
    }
    else {
        circos.genomicTrackPlotRegion(df, ylim = c(min(max(y_range)), max(y_range)),
            panel.fun = function(region, value, ...) {
                i = getI(...)
                circos.genomicLines(region, value, col = col[i],
                  lwd = lwd[i], lty = lty[i], type = type[i],
                  border = border[i], area = area[i], baseline = baseline[i])
            }, ...)
    }
}
environment(XY_circos.genomicDensity) <- asNamespace('circlize')

pseudo_bulk_seurat_mean_signac <- function(seurat_obj=seurat_obj,num_split=num_split,seed.use=seed.use,prefix=prefix,slot=slot,assay=assay){
  set.seed(seed.use)
  require(Seurat)
  genes.use <- rownames(seurat_obj)
  cell.sets1 <- split(colnames(seurat_obj), sort(1:length(colnames(seurat_obj))%%num_split))
  profile.set1 = matrix(, nrow = length(genes.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- GetAssayData(seurat_obj, slot = slot, assay = assay)[genes.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) mean(x)))
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix
    }
  }
  rownames(profile.set1) <- genes.use
  colnames(profile.set1) <- paste(prefix, 1:length(cell.sets1),sep="_")
  return(profile.set1)
}

pseudo_bulk_seurat_sum_signac <- function(seurat_obj=seurat_obj,num_split=num_split,seed.use=seed.use,prefix=prefix,slot=slot,assay=assay){
  set.seed(seed.use)
  require(Seurat)
  genes.use <- rownames(seurat_obj)
  cell.sets1 <- split(colnames(seurat_obj), sort(1:length(colnames(seurat_obj))%%num_split))
  profile.set1 = matrix(, nrow = length(genes.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- GetAssayData(seurat_obj, slot = slot, assay = assay)[genes.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix
    }
  }
  rownames(profile.set1) <- genes.use
  colnames(profile.set1) <- paste(prefix, 1:length(cell.sets1),sep="_")
  return(profile.set1)
}

XY_getGenePositions <- function (gene_names, ensembl_version = "dec2016.archive.ensembl.org",
    species = "human")
{
    if (species == "human") {
        ensembl = biomaRt::useMart("ensembl",dataset = "hsapiens_gene_ensembl")
        gene_positions <- biomaRt::getBM(attributes = c("ensembl_gene_id",
            "hgnc_symbol", "chromosome_name", "start_position",
            "end_position"), filters = "hgnc_symbol", values = gene_names,
            mart = ensembl)
    }
    else {
        ensembl = biomaRt::useMart("ensembl",dataset = "mmusculus_gene_ensembl")
        gene_positions <- biomaRt::getBM(attributes = c("ensembl_gene_id",
            "mgi_symbol", "chromosome_name", "start_position",
            "end_position"), filters = "mgi_symbol", values = gene_names,
            mart = ensembl)
    }
    gene_positions = gene_positions[!duplicated(gene_positions[,
        2]), ]
    gene_positions[which(gene_positions[, 3] == "X"), 3] = 23
    gene_positions[which(gene_positions[, 3] == "Y"), 3] = 24
    gene_positions[which(gene_positions[, 3] == "MT"), 3] = 0
    gene_positions[which(nchar(gene_positions[, 3]) > 2), 3] = 0
    gene_positions = gene_positions[order(as.numeric(gene_positions[,
        3]), decreasing = F), ]
    return(gene_positions)
}
environment(XY_getGenePositions) <- asNamespace('biomaRt')

XY_gene_atomization <- function (m) {
    tempzz = unlist(lapply(m@Abstract, function(x) {
        tempa = strsplit(x, ".  ", fixed = T)
        # tempa1 = which(nchar(tempa[[1]]) == max(nchar(tempa[[1]])))
        # tempb = unlist(strsplit(tempa[[1]][tempa1], ".", fixed = T))
        tempb = unlist(strsplit(unlist(tempa), ".", fixed = T))
        tempc = unlist(strsplit(tempb, ",", fixed = T))
        tempd = unlist(strsplit(tempc, ":", fixed = T))
        tempe = unlist(strsplit(tempd, ";", fixed = T))
        tempe1 = unlist(strsplit(tempe, "'", fixed = T))
        tempf = unlist(strsplit(tempe1, " ", fixed = T))
        tempf = unlist(strsplit(tempf, "-", fixed = T))
        return(tempf)
    }))
    tempi = as.data.frame(table(tempzz))
    tempj = unlist(lapply(common_words_new, function(x) {
        tempoo = which(as.character(tempi[, 1]) == x)
        if (length(tempoo) != 0)
            return(tempoo)
    }))
    tempk = tempi[-tempj, ]
    Class_genes <- as.character(HGNCdata$Approved.Symbol)
    Class_genes <- Class_genes[-which(Class_genes=="T")]
    Class_names <- as.character(HGNCdata$Approved.Name)
    Class_names <- Class_names[-which(Class_genes=="T")]
    templ = c(as.character(Class_genes))
    tempm = unlist(lapply(templ, function(x) {
        return(which(x == as.character(tempk$tempzz)))
    }))
    tempn = tempk[tempm, ]
    tempn$tempzz <- as.character(tempn$tempzz)
    tempn <- tempn[!duplicated(tempn$tempzz),]
    tempo = unlist(lapply(as.character(tempn$tempzz), function(x) {
        return(which(x == templ))
    }))
    Genes = as.character(Class_names[tempo])
    data_table = data.frame(Gene_symbol=as.character(tempn$tempzz),
      Genes=Genes,Freq=tempn$Freq)
    data_table$Gene_symbol <- as.character(data_table$Gene_symbol)
    data_table$Genes <- as.character(data_table$Genes)
    data_table = data_table[order(as.numeric(data_table$Freq), decreasing = T),]
    if (unique(data_table$Gene_symbol %in% "MLL3") ==TRUE) {
      data_table[data_table$Gene_symbol %in% "MLL3","Genes"] <- "mixed lineage leukemia 3"
    }
    write.table(data_table, file = "table.txt", sep = "\t", row.names = F)
    return(data_table)
  }
environment(XY_gene_atomization) <- asNamespace('pubmed.mineR')
